#pragma once
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <vector>

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

class NaiveStreamWriter : public IStreamWriter
{
public:
	NaiveStreamWriter(const std::filesystem::path& file_path);
	void write(const char *data, size_t size) override;
private:
	std::ofstream file;
};

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