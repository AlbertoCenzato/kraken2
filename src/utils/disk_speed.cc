#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <filesystem>

void measureMaxWriteSpeed(const std::filesystem::path& file_path) {
  std::vector<uint64_t> values = std::vector<uint64_t>(10 * 1024 * 1024);
  size_t total_written_bytes = 0;
  auto start = std::chrono::high_resolution_clock::now();
  std::ofstream output_file(file_path, std::ios::binary);
  for (uint64_t i = 0; i < 100; i++) {
    output_file.write(reinterpret_cast<const char *>(values.data()),
                      sizeof(uint64_t) * values.size());
    total_written_bytes += sizeof(uint64_t) * values.size();
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "Time to write file: "
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

void measureMaxReadSpeed(const std::filesystem::path& file_path) {
  std::vector<uint64_t> values = std::vector<uint64_t>(10 * 1024 * 1024);
  size_t total_read_bytes = 0;
  auto start = std::chrono::high_resolution_clock::now();
  std::ifstream input_file(file_path, std::ios::binary);
  while (input_file.read(reinterpret_cast<char *>(values.data()),
                         sizeof(uint64_t) * values.size())) {
    total_read_bytes += sizeof(uint64_t) * values.size();
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::cout << "Time to read file: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     start)
                   .count()
            << "ms" << std::endl;
  std::cout << "Read speed: "
            << (total_read_bytes / 1024 / 1024) /
                   (std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                          start)
                        .count() /
                    1000.0)
            << "MB/s" << std::endl;
}

int main() {
	auto filename = "output_file.bin";
	measureMaxWriteSpeed(filename);
	measureMaxReadSpeed(filename);
	std::filesystem::remove(filename);
	return 0;
}