#include "ListFiles_x64.h"

/*
@author:CodingMengmeng
@theme:获取指定文件夹下的所有文件名
@time:2017-1-13 11:46:22
@blog:http://www.cnblogs.com/codingmengmeng/
*/


void GetFilesWithType(string filePath, string subname, vector<string>& ofiles)
{
	vector<string> files;

	// 获取该路径下的所有文件  
	getFiles(filePath, files);

	string str1;
	for (int i = 0; i < files.size(); i++)
	{
		str1 = files[i].c_str();
		if (str1.find(subname) != string::npos)
		{
			ofiles.push_back(files[i]);
		}
	}


}


void getFiles(string path, vector<string>& files)
{
	//文件句柄  
	_int64 hFile = 0;
	//文件信息，声明一个存储文件信息的结构体  
	struct _finddatai64_t fileinfo;
	string p;//字符串，存放路径
	if ((hFile = _findfirsti64(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)//若查找成功，则进入
	{
		do
		{
			//如果是目录,迭代之（即文件夹内还有文件夹）  
			if ((fileinfo.attrib &  _A_SUBDIR))
			{
				continue;

				//文件名不等于"."&&文件名不等于".."
				//.表示当前目录
				//..表示当前目录的父目录
				//判断时，两者都要忽略，不然就无限递归跳不出去了！
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
				{
					getFiles(p.assign(path).append("\\").append(fileinfo.name), files);

				}
			}
			//如果不是,加入列表  
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
			}
		} while (_findnexti64(hFile, &fileinfo) == 0);
		//_findclose函数结束查找
		_findclose(hFile);
	}
}


/*

int main(){
char * filePath = "E:\\experiment data\\20171115 quality control 1qc\\with qc\\group1\\MMStack_Pos0";//自己设置目录
vector<string> files;

////获取该路径下的所有文件
getFiles(filePath, files);

char str[30];
int size = files.size();
for (int i = 0; i < size; i++)
{
cout << files[i].c_str() << endl;
}
}

*/