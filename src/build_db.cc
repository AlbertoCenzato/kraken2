/*
 * Copyright 2013-2023, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#include "compact_hash.h"
#include "kraken2_data.h"
#include "kraken2_headers.h"
#include "kv_store.h"
#include "merge_sort_files.h"
#include "mmscanner.h"
#include "seqreader.h"
#include "taxonomy.h"
#include "utilities.h"

#include <algorithm>
#include <filesystem>

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::map;
using std::ofstream;
using std::queue;
using std::string;
using std::vector;
using namespace kraken2;

// These will likely be overridden by wrapper script,
// just keeping sane defaults in case they aren't
#define DEFAULT_BLOCK_SIZE (10 * 1024 * 1024) // 10 MB
#define DEFAULT_SUBBLOCK_SIZE (1024)

struct Options {
  string ID_to_taxon_map_filename;
  string ncbi_taxonomy_directory;
  string hashtable_filename;
  string options_filename;
  string taxonomy_filename;
  size_t block_size;
  size_t subblock_size;
  size_t requested_bits_for_taxid;
  int num_threads;
  bool input_is_protein;
  ssize_t k, l;
  size_t capacity;
  size_t maximum_capacity;
  uint64_t spaced_seed_mask;
  uint64_t toggle_mask;
  uint64_t min_clear_hash_value;
  bool deterministic_build;
  string input_filename;
};

void ParseCommandLine(int argc, char **argv, Options &opts);
void usage(int exit_code = EX_USAGE);
vector<string> ExtractNCBISequenceIDs(const string &header);
void ProcessSequenceFast(const string &seq, taxid_t taxid,
                         CompactHashTable &hash, const Taxonomy &tax,
                         MinimizerScanner &scanner,
                         uint64_t min_clear_hash_value);
void ProcessSequence(const Options &opts, const string &seq, taxid_t taxid,
                     CompactHashTable &hash, const Taxonomy &tax);

std::vector<IdxValue> ProcessSequence2(const Options &opts, const string &seq,
                                       taxid_t taxid, const Taxonomy &tax,
                                       size_t value_bits, size_t capacity);

void ProcessSequencesFast(const Options &opts,
                          const map<string, taxid_t> &ID_to_taxon_map,
                          CompactHashTable &kraken_index,
                          const Taxonomy &taxonomy);
void ProcessSequences(const Options &opts,
                      const map<string, taxid_t> &ID_to_taxon_map,
                      CompactHashTable &kraken_index, const Taxonomy &taxonomy,
                      std::istream &input_stream, size_t value_bits,
                      size_t capacity);
void SetMinimizerLCA(CompactHashTable &hash, uint64_t minimizer, taxid_t taxid,
                     const Taxonomy &tax);
void ReadIDToTaxonMap(map<string, taxid_t> &id_map, string &filename);
void GenerateTaxonomy(Options &opts, map<string, taxid_t> &id_map);

struct TaxonSeqPair {
  taxid_t taxon;
  string seq;
};

int main(int argc, char **argv) {
  Options opts;
  opts.spaced_seed_mask = DEFAULT_SPACED_SEED_MASK;
  opts.toggle_mask = DEFAULT_TOGGLE_MASK;
  opts.input_is_protein = false;
  opts.num_threads = 1;
  opts.block_size = DEFAULT_BLOCK_SIZE;
  opts.subblock_size = DEFAULT_SUBBLOCK_SIZE;
  opts.requested_bits_for_taxid = 0;
  opts.min_clear_hash_value = 0;
  opts.maximum_capacity = 0;
  opts.deterministic_build = true;
  opts.input_filename =
      "/workspaces/kraken2/db/test/library/archaea/library.fna";
  ParseCommandLine(argc, argv, opts);

  omp_set_num_threads(opts.num_threads);

  map<string, taxid_t> ID_to_taxon_map;

  ReadIDToTaxonMap(ID_to_taxon_map, opts.ID_to_taxon_map_filename);
  GenerateTaxonomy(opts, ID_to_taxon_map);

  std::cerr << "Taxonomy parsed and converted2." << std::endl;

  Taxonomy taxonomy(opts.taxonomy_filename.c_str());
  taxonomy.GenerateExternalToInternalIDMap();
  size_t bits_needed_for_value = 1;
  while ((1 << bits_needed_for_value) < (ssize_t)taxonomy.node_count())
    bits_needed_for_value++;
  if (opts.requested_bits_for_taxid > 0 &&
      bits_needed_for_value > opts.requested_bits_for_taxid)
    errx(EX_DATAERR, "more bits required for storing taxid");

  size_t bits_for_taxid = bits_needed_for_value;
  if (bits_for_taxid < opts.requested_bits_for_taxid)
    bits_for_taxid = opts.requested_bits_for_taxid;

  auto actual_capacity = opts.capacity;
  if (opts.maximum_capacity) {
    double frac = opts.maximum_capacity * 1.0 / opts.capacity;
    if (frac > 1)
      errx(EX_DATAERR, "maximum capacity larger than requested capacity");
    opts.min_clear_hash_value = (uint64_t)((1 - frac) * UINT64_MAX);
    actual_capacity = opts.maximum_capacity;
  }

  CompactHashTable kraken_index(actual_capacity, 32 - bits_for_taxid,
                                bits_for_taxid);
  std::cerr << "CHT created with " << bits_for_taxid
            << " bits reserved for taxid." << std::endl;

  if (opts.deterministic_build) {
    if (opts.input_filename.empty()) {
      ProcessSequences(opts, ID_to_taxon_map, kraken_index, taxonomy, std::cin,
                       bits_for_taxid, actual_capacity);
    } else {
      ifstream input_file(opts.input_filename);
      if (!input_file.good())
        err(EX_NOINPUT, "unable to read from '%s'",
            opts.input_filename.c_str());
      ProcessSequences(opts, ID_to_taxon_map, kraken_index, taxonomy,
                       input_file, bits_for_taxid, actual_capacity);
    }
  } else
    ProcessSequencesFast(opts, ID_to_taxon_map, kraken_index, taxonomy);

  std::cerr << "Writing data to disk... " << std::flush;
  kraken_index.WriteTable(opts.hashtable_filename.c_str());

  IndexOptions index_opts;
  index_opts.k = opts.k;
  index_opts.l = opts.l;
  index_opts.spaced_seed_mask = opts.spaced_seed_mask;
  index_opts.toggle_mask = opts.toggle_mask;
  index_opts.dna_db = !opts.input_is_protein;
  index_opts.minimum_acceptable_hash_value = opts.min_clear_hash_value;
  index_opts.revcom_version = CURRENT_REVCOM_VERSION;
  index_opts.db_version = 0;
  index_opts.db_type = 0;
  ofstream opts_fs(opts.options_filename);
  opts_fs.write((char *)&index_opts, sizeof(index_opts));
  if (!opts_fs.good())
    errx(EX_OSERR, "Unable to write options file %s",
         opts.options_filename.c_str());
  opts_fs.close();
  std::cerr << " complete." << std::endl;

  return 0;
}

// A quick but nondeterministic build
void ProcessSequencesFast(const Options &opts,
                          const map<string, taxid_t> &ID_to_taxon_map,
                          CompactHashTable &kraken_index,
                          const Taxonomy &taxonomy) {
  size_t processed_seq_ct = 0;
  size_t processed_ch_ct = 0;

#pragma omp parallel
  {
    Sequence sequence;
    MinimizerScanner scanner(opts.k, opts.l, opts.spaced_seed_mask,
                             !opts.input_is_protein, opts.toggle_mask);

    BatchSequenceReader reader;

    while (true) {
      // Declaration of "ok" and break need to be done outside of critical
      // section to conform with OpenMP spec.
      bool ok;
#pragma omp critical(reader)
      ok = reader.LoadBlock(std::cin, opts.block_size);
      if (!ok)
        break;
      while (reader.NextSequence(sequence)) {
        auto all_sequence_ids = ExtractNCBISequenceIDs(sequence.header);
        taxid_t taxid = 0;
        for (auto &seqid : all_sequence_ids) {
          if (ID_to_taxon_map.count(seqid) == 0 ||
              ID_to_taxon_map.at(seqid) == 0)
            continue;
          auto ext_taxid = ID_to_taxon_map.at(seqid);
          taxid = taxonomy.LowestCommonAncestor(
              taxid, taxonomy.GetInternalID(ext_taxid));
        }
        if (taxid) {
          // Add terminator for protein sequences if not already there
          if (opts.input_is_protein && sequence.seq.back() != '*')
            sequence.seq.push_back('*');
          ProcessSequenceFast(sequence.seq, taxid, kraken_index, taxonomy,
                              scanner, opts.min_clear_hash_value);
#pragma omp atomic
          processed_seq_ct++;
#pragma omp atomic
          processed_ch_ct += sequence.seq.size();
        }
      }
      if (isatty(fileno(stderr))) {
#pragma omp critical(status_update)
        std::cerr << "\rProcessed " << processed_seq_ct << " sequences ("
                  << processed_ch_ct << " "
                  << (opts.input_is_protein ? "aa" : "bp") << ")...";
      }
    }
  }
  if (isatty(fileno(stderr)))
    std::cerr << "\r";
  std::cerr << "Completed processing of " << processed_seq_ct << " sequences, "
            << processed_ch_ct << " " << (opts.input_is_protein ? "aa" : "bp")
            << std::endl;
}

bool FileExists(const std::string &filename) {
  struct stat buffer;
  return (stat(filename.c_str(), &buffer) == 0);
}

// Slightly slower but deterministic when multithreaded
void ProcessSequences(const Options &opts,
                      const map<string, taxid_t> &ID_to_taxon_map,
                      CompactHashTable &kraken_index, const Taxonomy &taxonomy,
                      std::istream &input_stream, size_t value_bits,
                      size_t capacity) {
  auto start = std::chrono::steady_clock::now();
  size_t processed_seq_ct = 0;
  size_t processed_ch_ct = 0;

  Sequence sequence;
  BatchSequenceReader reader;

  std::filesystem::path output_hash_dir = "/workspaces/kraken2/db/test/hashes";
  std::filesystem::path output_merge_dir =
      "/workspaces/kraken2/db/test/merged_hashes";

  std::filesystem::create_directories(output_hash_dir);
  std::filesystem::create_directories(output_merge_dir);

  auto start_hashing = std::chrono::steady_clock::now();
  while (reader.LoadBlock(input_stream, DEFAULT_BLOCK_SIZE)) {
    std::cout << "Doing things..." << std::endl;
    while (reader.NextSequence(sequence)) {
      auto all_sequence_ids = ExtractNCBISequenceIDs(sequence.header);
      taxid_t taxid = 0;
      for (auto &seqid : all_sequence_ids) {
        if (ID_to_taxon_map.count(seqid) == 0 || ID_to_taxon_map.at(seqid) == 0)
          continue;
        auto ext_taxid = ID_to_taxon_map.at(seqid);
        taxid = taxonomy.LowestCommonAncestor(
            taxid, taxonomy.GetInternalID(ext_taxid));
      }
      if (taxid) {
        // Add terminator for protein sequences if not already there
        if (opts.input_is_protein && sequence.seq.back() != '*')
          sequence.seq.push_back('*');
        auto hashes = ProcessSequence2(opts, sequence.seq, taxid, taxonomy,
                                       value_bits, capacity);

        std::sort(hashes.begin(), hashes.end());

        uint64_t index = 0;
        auto output_path = output_hash_dir / (std::to_string(taxid) + "_" +
                                              std::to_string(index) + ".bin");
        while (FileExists(output_path)) {
          index++;
          output_path = output_hash_dir / (std::to_string(taxid) + "_" +
                                           std::to_string(index) + ".bin");
        }

        auto output_file = std::ofstream(output_path, std::ios::binary);
        output_file.write(reinterpret_cast<const char *>(hashes.data()),
                          hashes.size() * sizeof(IdxValue));

        // ProcessSequence(opts, sequence.seq, taxid, kraken_index, taxonomy);
        processed_seq_ct++;
        processed_ch_ct += sequence.seq.size();
      }
    }
    if (isatty(fileno(stderr))) {
      std::cerr << "\rProcessed " << processed_seq_ct << " sequences ("
                << processed_ch_ct << " "
                << (opts.input_is_protein ? "aa" : "bp") << ")...";
    }
  }
  auto end_hashing = std::chrono::steady_clock::now();
  std::cerr << "Hashing took "
            << std::chrono::duration_cast<std::chrono::seconds>(end_hashing -
                                                                start_hashing)
                   .count()
            << " seconds" << std::endl;

  auto start_merge = std::chrono::steady_clock::now();
  multiStepMerge(output_hash_dir, output_merge_dir,
                 output_merge_dir / "merge.bin");
  auto end_merge = std::chrono::steady_clock::now();
  std::cerr << "Merging took "
            << std::chrono::duration_cast<std::chrono::seconds>(end_merge -
                                                                start_merge)
                   .count()
            << " seconds" << std::endl;

  std::cout << "Value bits: " << value_bits << std::endl;
  std::cout << "Capacity: " << capacity << std::endl;

  auto start_hash_table = std::chrono::steady_clock::now();
  makeHashTable(output_merge_dir / "merge.bin",
                output_merge_dir / "hash_table.bin", capacity, 32 - value_bits,
                value_bits);
  auto end_hash_table = std::chrono::steady_clock::now();
  std::cerr << "Hash table creation took "
            << std::chrono::duration_cast<std::chrono::seconds>(
                   end_hash_table - start_hash_table)
                   .count()
            << " seconds" << std::endl;

  if (isatty(fileno(stderr)))
    std::cerr << "\r";
  std::cerr << "Completed processing of " << processed_seq_ct << " sequences, "
            << processed_ch_ct << " " << (opts.input_is_protein ? "aa" : "bp")
            << std::endl;
  auto end = std::chrono::steady_clock::now();
  std::cerr
      << "Processing took "
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << " seconds" << std::endl;
}

// This function exists to deal with NCBI's use of \x01 characters to denote
// the start of a new FASTA header in the same line (for non-redundant DBs).
// We return all sequence IDs in a header line, not just the first.
vector<string> ExtractNCBISequenceIDs(const string &header) {
  vector<string> list;
  string current_str;

  bool in_id = true;
  // start loop at first char after '>'
  for (size_t i = 1; i < header.size(); ++i) {
    if (header[i] == 0x01) {
      // 0x01 starts new ID at next char
      if (!current_str.empty())
        list.push_back(current_str);
      current_str.clear();
      in_id = true;
    } else if (in_id && isspace(header[i])) {
      // spaces end ID
      if (!current_str.empty())
        list.push_back(current_str);
      current_str.clear();
      in_id = false;
    } else if (in_id) {
      // Build ID string char by char
      current_str.push_back(header[i]);
    }
  }
  if (!current_str.empty())
    list.push_back(current_str);
  return list;
}

void SetMinimizerLCA(CompactHashTable &hash, uint64_t minimizer, taxid_t taxid,
                     const Taxonomy &tax) {
  hvalue_t old_value = 0;
  hvalue_t new_value = taxid;
  while (!hash.CompareAndSet(minimizer, new_value, &old_value))
    new_value = tax.LowestCommonAncestor(old_value, taxid);
}

void ProcessSequenceFast(const string &seq, taxid_t taxid,
                         CompactHashTable &hash, const Taxonomy &tax,
                         MinimizerScanner &scanner,
                         uint64_t min_clear_hash_value) {
  scanner.LoadSequence(seq);
  uint64_t *minimizer_ptr;
  while ((minimizer_ptr = scanner.NextMinimizer())) {
    if (scanner.is_ambiguous())
      continue;
    if (min_clear_hash_value &&
        MurmurHash3(*minimizer_ptr) < min_clear_hash_value)
      continue;
    hvalue_t existing_taxid = 0;
    hvalue_t new_taxid = taxid;
    while (!hash.CompareAndSet(*minimizer_ptr, new_taxid, &existing_taxid)) {
      new_taxid = tax.LowestCommonAncestor(new_taxid, existing_taxid);
    }
  }
}

IdxValue getIdxValue(hkey_t key, hvalue_t val, size_t value_bits,
                     size_t capacity) {
  uint64_t hc = MurmurHash3(key);
  size_t idx = hc % capacity;

  hkey_t compacted_key = hc >> (32 + value_bits);
  uint32_t value = computeCompactHash(compacted_key, val, value_bits);
  return IdxValue{.idx = idx, .value = value};
}

std::vector<IdxValue> ProcessSequence2(const Options &opts, const string &seq,
                                       taxid_t taxid, const Taxonomy &tax,
                                       size_t value_bits, size_t capacity) {
  const int set_ct = 256;
  // for each block in the sequence
  auto hashes = std::vector<IdxValue>{};
  for (size_t j = 0; j < seq.size(); j += opts.block_size) {
    // block: fixed length subsequence of DNA/protein, handled in series,
    //   consecutive blocks overlap by k-1 characters
    size_t block_start = j;
    size_t block_finish = j + opts.block_size + opts.k - 1;
    if (block_finish > seq.size())
      block_finish = seq.size();
    std::set<uint64_t> minimizer_sets[set_ct];

    // for each subblock in the block, gather minimizers in parallel
    for (size_t i = block_start; i < block_finish; i += opts.subblock_size) {
      // subblock: fixed length subsequence of block, handled in parallel,
      //   consecutive subblocks overlap by k-1 characters
      MinimizerScanner scanner(opts.k, opts.l, opts.spaced_seed_mask,
                               !opts.input_is_protein, opts.toggle_mask);
      size_t subblock_finish = i + opts.subblock_size + opts.k - 1;
      if (subblock_finish > block_finish)
        subblock_finish = block_finish;
      scanner.LoadSequence(seq, i, subblock_finish);
      uint64_t *minimizer_ptr;
      while ((minimizer_ptr = scanner.NextMinimizer())) {
        if (scanner.is_ambiguous())
          continue;
        auto hc = MurmurHash3(*minimizer_ptr);
        // Hash-based subsampling
        if (opts.min_clear_hash_value && hc < opts.min_clear_hash_value)
          continue;
        auto zone = hc % set_ct;
        minimizer_sets[zone].insert(*minimizer_ptr);
      }
    } // end subblock for loop

    // combine sets into sorted list
    std::vector<uint64_t> minimizer_lists[set_ct];
    size_t minimizer_list_prefix_sizes[set_ct + 1] = {0};
    for (int i = 0; i < set_ct; i++) {
      minimizer_lists[i].reserve(minimizer_sets[i].size());
      for (auto &m : minimizer_sets[i])
        minimizer_lists[i].push_back(m);
      sort(minimizer_lists[i].begin(), minimizer_lists[i].end());
      minimizer_list_prefix_sizes[i + 1] = minimizer_lists[i].size();
    }
    for (int i = 2; i <= set_ct; i++)
      minimizer_list_prefix_sizes[i] += minimizer_list_prefix_sizes[i - 1];
    std::vector<uint64_t> minimizer_list(minimizer_list_prefix_sizes[set_ct]);
    for (int i = 0; i < set_ct; i++) {
      std::copy(minimizer_lists[i].begin(), minimizer_lists[i].end(),
                minimizer_list.begin() + minimizer_list_prefix_sizes[i]);
    }

    for (auto minimizer : minimizer_list) {
      auto idx_val = getIdxValue(minimizer, taxid, value_bits, capacity);
      hashes.push_back(idx_val);
    }
  } // end block for loop
  return hashes;
}

void mergeSortStreams(std::vector<std::ifstream> &input_streams,
                      std::ostream &output_stream) {
  auto start = std::chrono::high_resolution_clock::now();
  size_t total_written_bytes = 0;
  std::vector<IdxValue> values(input_streams.size());
  std::vector<bool> eof(input_streams.size(), false);
  std::vector<bool> read_mask(input_streams.size(), true);
  while (!std::all_of(eof.begin(), eof.end(), [](bool x) { return x; })) {
    // read values from all files
    for (size_t i = 0; i < input_streams.size(); i++) {
      if (read_mask[i] && !eof[i]) {
        eof[i] = !input_streams[i].read(reinterpret_cast<char *>(&values[i]),
                                        sizeof(IdxValue));
        read_mask[i] = false;
      }
    }

    // find the smallest value
    const auto first_not_eof = std::find(eof.begin(), eof.end(), false);
    if (first_not_eof == eof.end()) {
      break;
    }

    size_t min_index = std::distance(eof.begin(), first_not_eof);
    for (size_t i = min_index + 1; i < input_streams.size(); i++) {
      if (!eof[i] && values[i] < values[min_index]) {
        min_index = i;
      }
    }

    read_mask[min_index] = true;
    const auto min_value = values[min_index];

    // write the smallest value to the output file
    output_stream.write(reinterpret_cast<const char *>(&min_value),
                        sizeof(IdxValue));
    total_written_bytes += sizeof(IdxValue);
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "Time to merge files: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << "ms" << std::endl;
  std::cout << "Write speed: "
            << (total_written_bytes / 1024 / 1024) /
                   (std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                          start)
                        .count() /
                    1000.0)
            << "MB/s" << std::endl;
}

void ProcessSequence(const Options &opts, const string &seq, taxid_t taxid,
                     CompactHashTable &hash, const Taxonomy &tax) {
  const int set_ct = 256;
  omp_lock_t locks[set_ct];
  for (int i = 0; i < set_ct; i++)
    omp_init_lock(&locks[i]);
  // for each block in the sequence
  for (size_t j = 0; j < seq.size(); j += opts.block_size) {
    // block: fixed length subsequence of DNA/protein, handled in series,
    //   consecutive blocks overlap by k-1 characters
    size_t block_start = j;
    size_t block_finish = j + opts.block_size + opts.k - 1;
    if (block_finish > seq.size())
      block_finish = seq.size();
    std::set<uint64_t> minimizer_sets[set_ct];

// for each subblock in the block, gather minimizers in parallel
#pragma omp parallel for schedule(dynamic)
    for (size_t i = block_start; i < block_finish; i += opts.subblock_size) {
      // subblock: fixed length subsequence of block, handled in parallel,
      //   consecutive subblocks overlap by k-1 characters
      MinimizerScanner scanner(opts.k, opts.l, opts.spaced_seed_mask,
                               !opts.input_is_protein, opts.toggle_mask);
      size_t subblock_finish = i + opts.subblock_size + opts.k - 1;
      if (subblock_finish > block_finish)
        subblock_finish = block_finish;
      scanner.LoadSequence(seq, i, subblock_finish);
      uint64_t *minimizer_ptr;
      while ((minimizer_ptr = scanner.NextMinimizer())) {
        if (scanner.is_ambiguous())
          continue;
        auto hc = MurmurHash3(*minimizer_ptr);
        // Hash-based subsampling
        if (opts.min_clear_hash_value && hc < opts.min_clear_hash_value)
          continue;
        auto zone = hc % set_ct;
        omp_set_lock(&locks[zone]);
        minimizer_sets[zone].insert(*minimizer_ptr);
        omp_unset_lock(&locks[zone]);
      }
    } // end subblock for loop

    // combine sets into sorted list
    std::vector<uint64_t> minimizer_lists[set_ct];
    size_t minimizer_list_prefix_sizes[set_ct + 1] = {0};
#pragma omp parallel for
    for (int i = 0; i < set_ct; i++) {
      minimizer_lists[i].reserve(minimizer_sets[i].size());
      for (auto &m : minimizer_sets[i])
        minimizer_lists[i].push_back(m);
      sort(minimizer_lists[i].begin(), minimizer_lists[i].end());
      minimizer_list_prefix_sizes[i + 1] = minimizer_lists[i].size();
    }
    for (int i = 2; i <= set_ct; i++)
      minimizer_list_prefix_sizes[i] += minimizer_list_prefix_sizes[i - 1];
    std::vector<uint64_t> minimizer_list(minimizer_list_prefix_sizes[set_ct]);
#pragma omp parallel for
    for (int i = 0; i < set_ct; i++) {
      std::copy(minimizer_lists[i].begin(), minimizer_lists[i].end(),
                minimizer_list.begin() + minimizer_list_prefix_sizes[i]);
    }

    size_t mm_ct = minimizer_list.size();
    std::vector<size_t> index_list(mm_ct);
    std::vector<bool> insertion_list(mm_ct);
    for (size_t i = 0; i < mm_ct; i++)
      insertion_list[i] = true;
    // Loop to enforce deterministic order
    while (!minimizer_list.empty()) {
// Gather insertion point information for all remaining minimizers
#pragma omp parallel for
      for (size_t i = 0; i < mm_ct; i++) {
        // once we've determined that a minimizer won't be an insertion,
        // don't bother calling FindIndex() again
        if (insertion_list[i]) {
          size_t idx;
          insertion_list[i] = !hash.FindIndex(minimizer_list[i], &idx);
          index_list[i] = idx;
        }
      }

      // Determine safe prefix of sorted set to insert in parallel
      std::set<uint64_t> novel_insertion_points;
      size_t safe_ct = 0;
      for (safe_ct = 0; safe_ct < mm_ct; safe_ct++) {
        if (insertion_list[safe_ct]) {
          if (novel_insertion_points.count(index_list[safe_ct]) > 0)
            break;
          novel_insertion_points.insert(index_list[safe_ct]);
        }
      }

// Adjust CHT values for all keys in the safe zone in parallel
#pragma omp parallel for
      for (size_t i = 0; i < safe_ct; i++) {
        SetMinimizerLCA(hash, minimizer_list[i], taxid, tax);
      }

      // Remove safe prefix and re-iterate to process remainder
      minimizer_list.erase(minimizer_list.begin(),
                           minimizer_list.begin() + safe_ct);
      index_list.erase(index_list.begin(), index_list.begin() + safe_ct);
      insertion_list.erase(insertion_list.begin(),
                           insertion_list.begin() + safe_ct);
      mm_ct -= safe_ct;
    } // end deterministic order while loop
  } // end block for loop

  for (int i = 0; i < set_ct; i++)
    omp_destroy_lock(&locks[i]);
}

void ReadIDToTaxonMap(map<string, taxid_t> &id_map, string &filename) {
  ifstream map_file(filename.c_str());
  if (!map_file.good())
    err(EX_NOINPUT, "unable to read from '%s'", filename.c_str());
  string line;

  while (getline(map_file, line)) {
    string seq_id;
    taxid_t taxid;
    istringstream iss(line);
    iss >> seq_id;
    iss >> taxid;
    if (taxid)
      id_map[seq_id] = taxid;
  }
}

void GenerateTaxonomy(Options &opts, map<string, taxid_t> &id_map) {
  NCBITaxonomy ncbi_taxonomy(opts.ncbi_taxonomy_directory + "/nodes.dmp",
                             opts.ncbi_taxonomy_directory + "/names.dmp");
  for (auto &kv_pair : id_map) {
    if (kv_pair.second != 0) {
      ncbi_taxonomy.MarkNode(kv_pair.second);
    }
  }
  ncbi_taxonomy.ConvertToKrakenTaxonomy(opts.taxonomy_filename.c_str());
}

void ParseCommandLine(int argc, char **argv, Options &opts) {
  int opt;
  long long sig;

  while ((opt = getopt(argc, argv, "?hB:b:c:FH:m:n:o:t:k:l:M:p:r:s:S:T:X")) !=
         -1) {
    switch (opt) {
    case 'h':
    case '?':
      usage(0);
      break;
    case 'B':
      sig = atoll(optarg);
      if (sig < 1)
        errx(EX_USAGE, "must have positive block size");
      opts.block_size = sig;
      break;
    case 'b':
      sig = atoll(optarg);
      if (sig < 1)
        errx(EX_USAGE, "must have positive subblock size");
      opts.subblock_size = sig;
      break;
    case 'r':
      sig = atoll(optarg);
      if (sig < 0)
        errx(EX_USAGE, "can't have negative bit storage");
      if (sig > 31)
        errx(EX_USAGE, "can't have more than 31 bits of storage for taxid");
      opts.requested_bits_for_taxid = sig;
      break;
    case 'p':
      opts.num_threads = atoll(optarg);
      if (opts.num_threads < 1)
        errx(EX_USAGE, "can't have negative number of threads");
      if (opts.num_threads > omp_get_max_threads())
        errx(EX_USAGE, "OMP only wants you to use %d threads",
             omp_get_max_threads());
      break;
    case 'H':
      opts.hashtable_filename = optarg;
      break;
    case 'm':
      opts.ID_to_taxon_map_filename = optarg;
      break;
    case 'n':
      opts.ncbi_taxonomy_directory = optarg;
      break;
    case 'o':
      opts.options_filename = optarg;
      break;
    case 't':
      opts.taxonomy_filename = optarg;
      break;
    case 'S':
      opts.spaced_seed_mask = strtol(optarg, nullptr, 2);
      break;
    case 'T':
      opts.toggle_mask = strtol(optarg, nullptr, 2);
      break;
    case 'k':
      sig = atoll(optarg);
      if (sig < 1)
        errx(EX_USAGE, "k must be positive integer");
      opts.k = sig;
      break;
    case 'l':
      sig = atoll(optarg);
      if (sig < 1)
        errx(EX_USAGE, "l must be positive integer");
      if (sig > 31)
        errx(EX_USAGE, "l must be no more than 31");
      opts.l = sig;
      break;
    case 'c':
      sig = atoll(optarg);
      if (sig < 1)
        errx(EX_USAGE, "capacity must be positive integer");
      opts.capacity = sig;
      break;
    case 'M':
      sig = atoll(optarg);
      if (sig < 1)
        errx(EX_USAGE, "max capacity must be positive integer");
      opts.maximum_capacity = sig;
      break;
    case 'F':
      opts.deterministic_build = false;
      break;
    case 'X':
      opts.input_is_protein = true;
      break;
    }
  }

  if (opts.spaced_seed_mask != DEFAULT_SPACED_SEED_MASK)
    ExpandSpacedSeedMask(opts.spaced_seed_mask, opts.input_is_protein
                                                    ? BITS_PER_CHAR_PRO
                                                    : BITS_PER_CHAR_DNA);
  if (opts.hashtable_filename.empty() ||
      opts.ID_to_taxon_map_filename.empty() ||
      opts.ncbi_taxonomy_directory.empty() || opts.options_filename.empty() ||
      opts.taxonomy_filename.empty()) {
    cerr << "missing mandatory filename parameter" << endl;
    usage();
  }
  if (opts.k == 0 || opts.l == 0 || opts.capacity == 0) {
    cerr << "missing mandatory integer parameter" << endl;
    usage();
  }
  if (opts.k < opts.l) {
    cerr << "k cannot be less than l" << endl;
    usage();
  }
  if (opts.block_size < opts.subblock_size) {
    cerr << "block size cannot be less than subblock size\n";
    usage();
  }
  if (opts.maximum_capacity > opts.capacity) {
    cerr << "maximum capacity option shouldn't specify larger capacity than "
            "normal"
         << endl;
    usage();
  }
}

void usage(int exit_code) {
  cerr << "Usage: build_db <options>\n"
       << "\n"
       << "Options (*mandatory):\n"
       << "* -H FILENAME   Kraken 2 hash table filename\n"
       << "* -m FILENAME   Sequence ID to taxon map filename\n"
       << "* -t FILENAME   Kraken 2 taxonomy filename\n"
       << "* -n DIR        NCBI taxonomy directory name\n"
       << "* -o FILENAME   Kraken 2 options filename\n"
       << "* -k INT        Set length of k-mers\n"
       << "* -l INT        Set length of minimizers\n"
       << "* -c INT        Set capacity of hash table\n"
       << "  -M INT        Set maximum capacity of hash table (MiniKraken)\n"
       << "  -S BITSTRING  Spaced seed mask\n"
       << "  -T BITSTRING  Minimizer toggle mask\n"
       << "  -X            Input seqs. are proteins\n"
       << "  -p INT        Number of threads\n"
       << "  -F            Use fast, nondeterministic building method\n"
       << "  -B INT        Read block size\n"
       << "  -b INT        Read subblock size\n"
       << "  -r INT        Bit storage requested for taxid" << endl;
  exit(exit_code);
}
