#pragma once
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <vector>
#include <array>
#include <mutex>
#include <thread>

class IStreamWriter
{
  public:
	virtual ~IStreamWriter() = default;
  virtual void write(const char *data, size_t size) = 0;
};

class IStreamReader
{
  public:
	virtual ~IStreamReader() = default;
  virtual size_t read(char *data, size_t size) = 0;
	virtual bool endOfStream() = 0;
};

class NaiveStreamReader : public IStreamReader
{
public:
	NaiveStreamReader(const std::filesystem::path& file_path);
	size_t read(char *data, size_t size) override;
	bool endOfStream() override;
private:
	std::ifstream file;
};

class BufferedStreamReader : public IStreamReader
{
public:
	BufferedStreamReader(const std::filesystem::path& file_path, size_t buffer_size);
	size_t read(char *data, size_t size) override;
	bool endOfStream() override;
private:
	std::vector<char> buffer;
	std::ifstream file;
};

class AsyncPrefetchStreamReader : public IStreamReader
{
public:
	AsyncPrefetchStreamReader(
		const std::filesystem::path& file_path, size_t buffer_size);
	size_t read(char *data, size_t size) override;
	bool endOfStream() override;
private:
	size_t buffer_size;
	std::mutex mutex;
	std::array<std::vector<char>,2> buffers;
	std::vector<char> *prefetch_buffer;
	std::vector<char> *active_read_buffer;
	size_t read_idx;
	std::ifstream file;
	std::jthread thread;
};

class NaiveStreamWriter : public IStreamWriter
{
public:
	NaiveStreamWriter(const std::filesystem::path& file_path);
	void write(const char *data, size_t size) override;
private:
	std::ofstream file;
};

class BufferedStreamWriter : public IStreamWriter
{
public:
	BufferedStreamWriter(const std::filesystem::path& file_path, size_t buffer_size);
	void write(const char *data, size_t size) override;
private:
	std::vector<char> buffer;
	std::ofstream file;
};

/*
class AsyncStreamWriter : public IStreamWriter
{
public:
	AsyncStreamWriter(const std::filesystem::path& file_path, size_t buffer_size);
	void write(const char *data, size_t size) override;
private:
	size_t buffer_size;
	std::mutex mutex;
	std::vector<char> buffer;
	std::ofstream file;
	std::jthread thread;
};
*/


class TestStreamReader : public IStreamReader
{
public:
	TestStreamReader(std::vector<char> data);
	size_t read(char *data, size_t size) override;
	bool endOfStream() override;
private:
	std::vector<char> data;
	size_t pos;
};

class TestStreamWriter : public IStreamWriter
{
public:
	std::vector<char> data;
	void write(const char *data, size_t size) override;
};

class NullWriter : public IStreamWriter
{
public:
	void write(const char *data, size_t size) override;
};