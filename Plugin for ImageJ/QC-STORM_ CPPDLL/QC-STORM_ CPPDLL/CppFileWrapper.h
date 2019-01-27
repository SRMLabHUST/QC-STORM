/*
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU LESSER GENERAL PUBLIC LICENSE for more details.

You should have received a copy of the GNU LESSER GENERAL PUBLIC LICENSE
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

#include <fstream>      // std::ifstream, std::ofstream

#include <iostream>
#include <string>
using namespace std;




class CppBinaryReadFile {
private:
	string FileName;
	std::ifstream infile;

	long long FileLengh;

public:

	CppBinaryReadFile(string filename)
	{
		this->FileName = filename;
		infile.open(filename, std::ifstream::binary);


		// get length of file:
		infile.seekg(0, ios::end);
		FileLengh = infile.tellg();
		infile.seekg(0, ios::beg);
	}

	void Read(void *buffer, long long length)
	{
		infile.read((char*)buffer, length);

	}

	bool IsOpen()
	{
		return infile.is_open();
	}

	//ios::beg, ios::cur, ios::end
	void Seek(int pos, ios_base::seekdir SeekMode)
	{
		infile.seekg(pos, SeekMode);
	}

	long long Tell()
	{
		return infile.tellg();
	}

	long long GetLength()
	{
		return FileLengh;
	}

	void Close()
	{
		infile.close();
	}
};


class CppBinaryWriteFile {

private:
	string FileName;
	std::ofstream outfile;

public:

	CppBinaryWriteFile(string filename)
	{
		this->FileName = filename;

		outfile.open(filename, std::ofstream::binary);
	}

	//ios::beg, ios::cur, ios::end
	void Seek(int pos, ios_base::seekdir SeekMode)
	{
		outfile.seekp(pos, SeekMode);
	}

	long long Tell()
	{
		return outfile.tellp();
	}

	void Write(void* buffer, long long length)
	{
		// write to outfile
		outfile.write((char*)buffer, length);

	}

	void Close()
	{
		outfile.close();

	}

};



#include <iostream>       // std::cout, std::hex
#include <string>         // std::string, std::u32string
#include <locale>         // std::wstring_convert
#include <codecvt>        // std::codecvt_utf8
#include <cstdint>        // std::uint_least32_t


class EncodeConvert {
public:

	static wstring s2ws(const std::string& str)
	{
		using convert_typeX = std::codecvt_utf8<wchar_t>;
		std::wstring_convert<convert_typeX, wchar_t> converterX;

		return converterX.from_bytes(str);
	}

	static string ws2s(const std::wstring& wstr)
	{
		using convert_typeX = std::codecvt_utf8<wchar_t>;
		std::wstring_convert<convert_typeX, wchar_t> converterX;

		return converterX.to_bytes(wstr);
	}

};


