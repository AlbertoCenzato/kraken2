#include "streams.h"
#include <iostream>

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

BufferedStreamReader::BufferedStreamReader(
	const std::filesystem::path& file_path, size_t buffer_size)
	: buffer(buffer_size)
{
	file.rdbuf()->pubsetbuf(buffer.data(), buffer.size());
	file.open(file_path, std::ios::binary);
}

size_t BufferedStreamReader::read(char *data, size_t size) 
{
	file.read(data, size);
	return file.gcount();
}

bool BufferedStreamReader::endOfStream() 
{
	return file.eof();
}

AsyncPrefetchStreamReader::AsyncPrefetchStreamReader(
		const std::filesystem::path& file_path, size_t buffer_size)
	: buffer_size{buffer_size}, mutex{}, 
		buffers{std::vector<char>(), std::vector<char>()},
		prefetch_buffer{&buffers[0]}, active_read_buffer{&buffers[1]}, read_idx{0}
{
	file.rdbuf()->pubsetbuf(prefetch_buffer->data(), prefetch_buffer->size());
	file.open(file_path, std::ios::binary);
	thread = std::jthread([this](std::stop_token st) {
		while (!st.stop_requested()) {
			auto lock = std::lock_guard(this->mutex);
			if (!this->prefetch_buffer->empty()) {
				// NOTE(alberto): buffer is not empty, wait for it to be consumed by AsyncPrefetchStreamReader::read()
				continue;
			}

			if (file.eof()) {
				// NOTE(alberto): should prefetch some data but end of file was reached
				// so we set the buffer to nullptr to signal the reader that it's done
				this->prefetch_buffer = nullptr;
				break;
			}

			this->prefetch_buffer->resize(this->buffer_size);
			file.read(this->prefetch_buffer->data(), this->prefetch_buffer->size());
			//std::cout << "Prefetched " << file.gcount() << " bytes" << std::endl;
		}
	});
}

size_t AsyncPrefetchStreamReader::read(char *data, size_t size) 
{
	while (read_idx >= active_read_buffer->size()) {
		read_idx = 0;
		active_read_buffer->clear();
		{
			auto lock = std::lock_guard(this->mutex);
			std::swap(prefetch_buffer, active_read_buffer);
		}

		if (endOfStream()) {
			return 0;
		}
	}

	size_t to_read = std::min(size, active_read_buffer->size() - read_idx);
	const auto start = active_read_buffer->begin() + read_idx;
	std::copy(start, start + to_read, data);
	read_idx += to_read;
	return to_read;
}

bool AsyncPrefetchStreamReader::endOfStream()
{
	return !active_read_buffer;
}

NaiveStreamWriter::NaiveStreamWriter(const std::filesystem::path& file_path) 
	: file{file_path, std::ios::binary} 
{

}

void NaiveStreamWriter::write(const char *data, size_t size) 
{
	file.write(data, size);
}

BufferedStreamWriter::BufferedStreamWriter(
	const std::filesystem::path& file_path, size_t buffer_size)
	: buffer(buffer_size)
{
	file.rdbuf()->pubsetbuf(buffer.data(), buffer.size());
	file.open(file_path, std::ios::binary);
}

void BufferedStreamWriter::write(const char *data, size_t size)
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

void NullWriter::write(const char *data, size_t size) 
{
	// do nothing
}
