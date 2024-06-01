#include "streams.h"

NaiveStreamReader::NaiveStreamReader(const std::filesystem::path& file_path) 
	: file{file_path, std::ios::binary} 
{

}

size_t NaiveStreamReader::read(char *data, size_t size)
{
	file.read(data, size);
	return file.gcount();
}

bool NaiveStreamReader::endOfStream()
{
	return file.eof();
}

NaiveStreamWriter::NaiveStreamWriter(const std::filesystem::path& file_path) 
	: file{file_path, std::ios::binary} 
{

}

void NaiveStreamWriter::write(const char *data, size_t size) 
{
	file.write(data, size);
}


TestStreamReader::TestStreamReader(std::vector<char> data) 
	: data{std::move(data)}, pos{0}
{

}
	
size_t TestStreamReader::read(char *data, size_t size) 
{
	size_t to_read = std::min(size, this->data.size() - pos);
	std::copy(this->data.begin() + pos, this->data.begin() + pos + to_read, data);
	pos += to_read;
	return to_read;
}

bool TestStreamReader::endOfStream() 
{
	return pos >= data.size();
}

void TestStreamWriter::write(const char *data, size_t size) 
{
	this->data.insert(this->data.end(), data, data + size);
}