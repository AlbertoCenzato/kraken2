#include "merge_sort_files.h"
#include "utils/streams.h"
#include <iostream>

int main() {
  //std::filesystem::path src = "/workspaces/kraken2/db/test/hashes/76_0.bin";
  //std::filesystem::path dst = "/workspaces/kraken2/db/test/merged_hashes/76_0.bin";
	//makeHashTable(src, dst, 898584, 32, 32);

	multiStepMerge(
		"/mnt/data/acenzato/kraken2/db/test/hashes", 
		"/mnt/data/acenzato/kraken2/db/test/merged_hashes", 
		"/mnt/data/acenzato/kraken2/db/test/merged_hashes/merge.bin");

	return 0;
}