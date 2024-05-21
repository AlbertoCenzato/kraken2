#include "merge_sort_files.h"

int main() {
  std::filesystem::path src = "/workspaces/kraken2/db/test/hashes/76_0.bin";
  std::filesystem::path dst = "/workspaces/kraken2/db/test/merged_hashes/76_0.bin";

	makeHashTable(src, dst, 898584, 32, 32);

	return 0;
}