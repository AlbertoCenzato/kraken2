#include "streams.h"
#include <iostream>
#include <chrono>

int main() {
	auto path = "/mnt/data/acenzato/kraken2/db/test/merged_hashes/merge.bin";
	auto stream = BufferedStreamReader(path, 1024*1024);

	auto start = std::chrono::high_resolution_clock::now();
	size_t total = 0;
	while(!stream.endOfStream()) {
		//std::cout<< "Reading..." << std::endl;
		char buffer[1024];
		size_t bytes = stream.read(buffer, 1);
		//std::cout << "Read " << bytes << std::endl;
		total += bytes;
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Time to read file: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
	std::cout << "Read speed: " << (total / 1024 / 1024) / (std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0) << "MB/s" << std::endl;
}