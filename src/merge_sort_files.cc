#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <thread>
#include <utility>
#include <vector>

#include "merge_sort_files.h"

namespace fs = std::filesystem;

// order by idx and then by value
bool operator<(const IdxValue &lhs, const IdxValue &rhs) {
  return lhs.idx < rhs.idx || (lhs.idx == rhs.idx && lhs.value < rhs.value);
}

bool operator==(const IdxValue &lhs, const IdxValue &rhs) {
  return lhs.idx == rhs.idx && lhs.value == rhs.value;
}

struct Data {
  std::ifstream file_stream;
  IdxValue value;
  bool eof;
  bool read_mask;
};

void mergeSortStreams(std::vector<fs::path> &input_files,
                      std::ostream &output_stream) {
  auto start = std::chrono::high_resolution_clock::now();

  size_t total_processed_bytes = 0;
  std::vector<Data> data;
  for (const auto &path : input_files) {
    data.push_back(Data{.file_stream = std::ifstream{path, std::ios::binary},
                        .value = {},
                        .eof = false,
                        .read_mask = true});
  }

  auto last_value_written = IdxValue{0, 0};
  while (!std::all_of(data.begin(), data.end(),
                      [](const auto &x) { return x.eof; })) {
    // read values from all files
    for (size_t i = 0; i < data.size(); i++) {
      auto &d = data[i];
      if (d.read_mask && !d.eof) {
        d.file_stream.read(reinterpret_cast<char *>(&d.value),
                           sizeof(IdxValue));
        d.read_mask = false;
        d.eof = d.file_stream.eof();
      }
    }

    // find the smallest value
    const auto first_not_eof = std::find_if(
        data.begin(), data.end(), [](const auto &x) { return !x.eof; });
    if (first_not_eof == data.end()) {
      break;
    }

    auto min_iter = first_not_eof;
    for (auto iter = first_not_eof + 1; iter != data.end(); iter++) {
      if (!iter->eof && iter->value < min_iter->value) {
        min_iter = iter;
      }
    }

    min_iter->read_mask = true;
    const auto min_value = min_iter->value;

    total_processed_bytes += sizeof(IdxValue);
    if (min_value == last_value_written) {
      continue;
    }

    // write the smallest value to the output file
    output_stream.write(reinterpret_cast<const char *>(&min_value),
                        sizeof(IdxValue));

    last_value_written = min_value;
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "Time to merge files: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << "ms" << std::endl;
  std::cout << "Processing speed: "
            << (total_processed_bytes / 1024.0 / 1024.0) /
                   (std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                          start)
                        .count() /
                    1000.0)
            << "MB/s" << std::endl;
}

void multiStepMerge(const fs::path &hash_dir, const fs::path &merge_dir,
                    const fs::path &output_file) {
  auto start_merge = std::chrono::steady_clock::now();
  auto src_dir = hash_dir;
  fs::directory_iterator it(src_dir);
  auto file_paths =
      std::vector<fs::directory_entry>{it, fs::directory_iterator{}};
  size_t iteration = 0;
  while (file_paths.size() > 1) {
    std::cout << "Iteration " << iteration << "; files: " << file_paths.size()
              << std::endl;
    auto dst_dir = merge_dir / ("iteration_" + std::to_string(iteration));
    fs::create_directories(dst_dir);

    constexpr size_t batch_size = 64;
    for (size_t i = 0; i < file_paths.size(); i += batch_size) {
      std::cout << "Merged files " << i << std::endl;
      size_t files_in_batch = std::min(batch_size, file_paths.size() - i);
      auto batch_file_paths = std::vector<fs::path>();
      for (size_t j = 0; j < files_in_batch; j++) {
        batch_file_paths.push_back(file_paths[i + j].path());
      }

      auto output_path =
          dst_dir / ("merged_" + std::to_string(i / batch_size) + ".bin");
      auto output_stream = std::ofstream(output_path, std::ios::binary);
      mergeSortStreams(batch_file_paths, output_stream);
    }

    // remove old files
    if (iteration > 0) {
      auto prev_dst_dir =
          merge_dir / ("iteration_" + std::to_string(iteration - 1));
      fs::remove_all(prev_dst_dir);
    }
    src_dir = std::move(dst_dir);
    fs::directory_iterator it(src_dir);
    file_paths = std::vector<fs::directory_entry>{it, fs::directory_iterator{}};
    iteration++;
  }

  // move the last file to the output file
  if (file_paths.size() != 1) {
    std::cerr << "Error: expected 1 file, got " << file_paths.size()
              << std::endl;
    return;
  }
  fs::rename(file_paths[0].path(), output_file);

  auto end_merge = std::chrono::steady_clock::now();
  std::cerr << "Merging took "
            << std::chrono::duration_cast<std::chrono::seconds>(end_merge -
                                                                start_merge)
                   .count()
            << " seconds" << std::endl;
}

std::pair<std::vector<std::ifstream>, std::ofstream> createTestFiles() {
  for (int i = 0; i < 2; i++) {
    auto start = std::chrono::high_resolution_clock::now();
    std::ofstream input_file("input_file" + std::to_string(i) + ".bin",
                             std::ios::binary);
    uint64_t max = i == 0 ? size_t(1e8) : size_t(1e8) + 1;
    for (uint64_t j = 0; j < max; j++) {
      auto value = 2 * j + i;
      input_file.write(reinterpret_cast<const char *>(&value),
                       sizeof(uint64_t));
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Time to write file " << i << ": "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                       start)
                     .count()
              << "ms" << std::endl;
    std::cout << "Write speed: "
              << (max * sizeof(uint64_t) / 1024 / 1024) /
                     (std::chrono::duration_cast<std::chrono::milliseconds>(
                          end - start)
                          .count() /
                      1000.0)
              << "MB/s" << std::endl;
  }

  std::ofstream output_file("output_file.bin", std::ios::binary);
  auto files = std::vector<std::ifstream>();
  files.push_back(std::ifstream("input_file0.bin", std::ios::binary));
  files.push_back(std::ifstream("input_file1.bin", std::ios::binary));
  return std::make_pair(std::move(files), std::move(output_file));
}

template <typename T>
bool writeToFile(std::ofstream &file, std::streampos position, const T &data) {
  // Move the file pointer to the specified position
  file.seekp(position);
  if (!file) {
    std::cerr << "Failed to seek to position: " << position << std::endl;
    return false;
  }

  // Write the byte at the specified position
  file.write(reinterpret_cast<const char *>(&data), sizeof(T));
  if (!file) {
    std::cerr << "Failed to write byte to file" << std::endl;
    return false;
  }

  return true;
}

void makeHashTable(const std::filesystem::path &merged_hash_file,
                   const std::filesystem::path &output_file, size_t capacity,
                   size_t key_bits, size_t value_bits) {
  {
    auto hash_table_file = std::ofstream(output_file, std::ios::binary);
    hash_table_file.write((char *)&capacity, sizeof(capacity));

    // TODO(alberto): will be filled later
    size_t size = 0;
    hash_table_file.write((char *)&size, sizeof(size));
    hash_table_file.write((char *)&key_bits, sizeof(key_bits));
    hash_table_file.write((char *)&value_bits, sizeof(value_bits));
  }

  // NOTE(alberto): resize file to max capacity and fill with 0s
  const auto header_size =
      sizeof(capacity) + sizeof(size_t) + sizeof(key_bits) + sizeof(value_bits);
  const auto filesize = header_size + capacity * sizeof(IdxValue::value);
  std::filesystem::resize_file(output_file, filesize);

  auto input_file = std::ifstream(merged_hash_file, std::ios::binary);

  // NOTE(alberto): avoid destroying the file
  auto flags = std::ios::binary | std::ios::in | std::ios::out;
  auto hash_table_file = std::ofstream(output_file, flags);

  size_t size = 0;
  auto buffer = std::vector<IdxValue>(1024);
  while (true) {
    input_file.read((char *)buffer.data(), sizeof(IdxValue) * buffer.size());
    if (!input_file) {
      break;
    }

    for (const auto &idx_value : buffer) {
      if ((size % 1000000) == 0) {
        std::cout << "Processed " << size << " values ("
                  << size * sizeof(IdxValue) / (1024 * 1024) << "MB)"
                  << std::endl;
      }
      auto byte_idx = idx_value.idx * sizeof(idx_value.value) + header_size;
      if (byte_idx >= filesize) {
        std::cerr << "Error: byte index " << byte_idx << " is out of bounds"
                  << std::endl;
        return;
      }
      writeToFile(hash_table_file, byte_idx, idx_value.value);
      size++;
    }
  }

  // write the size of the hash table
  writeToFile(hash_table_file, sizeof(capacity), size);
}

/*
int main(int argc, char *argv[]) {
         if (argc < 2) {
                std::cerr << "Usage: " << argv[0] << " <input file> <input file>
[<input file>...] <output file>" << std::endl; return 1;
        }

        std::cout << "Merging " << argc - 2 << " files into " << argv[argc - 1]
<< std::endl; int num_input_files = argc - 2; std::ofstream
output_file(argv[argc - 1], std::ios::binary); if (!output_file.is_open()) {
                std::cerr << "Failed to open output file " << argv[argc - 1] <<
std::endl; return 1;
        }

        // open all files
        auto files = std::vector<std::ifstream>();
        for (int i = 0; i < num_input_files; i++) {
                files.push_back(std::ifstream(argv[i + 1], std::ios::binary));
                if (!files[i].is_open()) {
                        std::cerr << "Failed to open input file " << argv[i + 1]
<< std::endl; return 1;
                }
        }

        std::cout << "merging files" << std::endl;
        mergeSortStreams(files, output_file);

        fs::path output_hash_dir = "/workspaces/kraken2/db/test/hashes";
  fs::path output_merge_dir = "/workspaces/kraken2/db/test/merged_hashes";
        fs::path output_file =
"/workspaces/kraken2/db/test/merged_hashes/merged.bin";
        multiStepMerge(output_hash_dir, output_merge_dir, output_file);

        return 0;
}
*/
