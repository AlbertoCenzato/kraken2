#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>
#include <chrono>

bool lessThan(uint64_t a, uint64_t b) {
	return a < b;
}

void mergeSortStreams(std::vector<std::ifstream>& input_streams, std::ostream& output_stream)
{
	auto start = std::chrono::high_resolution_clock::now();
	size_t total_written_bytes = 0;
	std::vector<uint64_t> values(input_streams.size());
	std::vector<bool> eof(input_streams.size(), false);
	std::vector<bool> read_mask(input_streams.size(), true);
	while (!std::all_of(eof.begin(), eof.end(), [](bool x) { return x; })) {
		// read values from all files
		for (size_t i = 0; i < input_streams.size(); i++) {
			if (read_mask[i] && !eof[i]) {
				eof[i] = !input_streams[i].read(reinterpret_cast<char *>(&values[i]), sizeof(uint64_t));
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
			if (!eof[i] && lessThan(values[i], values[min_index])) {
				min_index = i;
			}
		}
		
		read_mask[min_index] = true;
		const auto min_value = values[min_index];

		// write the smallest value to the output file
		output_stream.write(reinterpret_cast<const char *>(&min_value), sizeof(uint64_t));
		total_written_bytes += sizeof(uint64_t);
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Time to merge files: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
	std::cout << "Write speed: " << (total_written_bytes / 1024 / 1024) / (std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0) << "MB/s" << std::endl;
}

void measureMaxWriteSpeed() {
	std::vector<uint64_t> values = std::vector<uint64_t>(1024*1024);
	size_t total_written_bytes = 0;
	auto start = std::chrono::high_resolution_clock::now();
	std::ofstream output_file("output_file.bin", std::ios::binary);
	for (uint64_t i = 0; i < 100; i++) {
		output_file.write(reinterpret_cast<const char *>(values.data()), sizeof(uint64_t)*values.size());
		total_written_bytes += sizeof(uint64_t)*values.size();
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Time to write file: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
	std::cout << "Write speed: " << (total_written_bytes / 1024 / 1024) / (std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0) << "MB/s" << std::endl;

}

std::pair<std::vector<std::ifstream>, std::ofstream> createTestFiles()
{
	for (int i = 0; i < 2; i++){
		auto start = std::chrono::high_resolution_clock::now();
		std::ofstream input_file("input_file" + std::to_string(i) + ".bin", std::ios::binary);
		uint64_t max = i == 0 ? size_t(1e8) : size_t(1e8) + 1;
		for (uint64_t j = 0; j < max; j++) {
			auto value = 2*j + i;
			input_file.write(reinterpret_cast<const char *>(&value), sizeof(uint64_t));
		}
		auto end = std::chrono::high_resolution_clock::now();
		std::cout << "Time to write file " << i << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
		std::cout << "Write speed: " << (max * sizeof(uint64_t) / 1024 / 1024) / (std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0) << "MB/s" << std::endl;
	}

	std::ofstream output_file("output_file.bin", std::ios::binary);
	auto files = std::vector<std::ifstream>();
	files.push_back(std::ifstream("input_file0.bin", std::ios::binary));
	files.push_back(std::ifstream("input_file1.bin", std::ios::binary));
	return std::make_pair(std::move(files), std::move(output_file));
}



int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <input file> <input file> [<input file>...] <output file>" << std::endl;
		return 1;
	}

	std::cout << "Merging " << argc - 2 << " files into " << argv[argc - 1] << std::endl;
	int num_input_files = argc - 2;
	std::ofstream output_file(argv[argc - 1], std::ios::binary);
	if (!output_file.is_open()) {
		std::cerr << "Failed to open output file " << argv[argc - 1] << std::endl;
		return 1;
	}

	// open all files
	auto files = std::vector<std::ifstream>();
	for (int i = 0; i < num_input_files; i++) {
		files.push_back(std::ifstream(argv[i + 1], std::ios::binary));
		if (!files[i].is_open()) {
			std::cerr << "Failed to open input file " << argv[i + 1] << std::endl;
			return 1;
		}
	}

	std::cout << "merging files" << std::endl;
	mergeSortStreams(files, output_file);

	return 0;
}