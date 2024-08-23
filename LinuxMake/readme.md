# Instructions for linux users:

1) Download the **Makefile** and create a new folder called **classes**. Therein, create two folders **src** and **include**.
   * According to user feedback you might need a third folder **obj**, please try this if compilation does not work.
   * For newer versions of C++: Add link flag -ldl in the Makefile.

2) Download all files from folder **Sources** and add them as follows to your folder structure:
   * main.cpp is on the same level as folder *classes*
   * *include* contains all header files, i.e. main.h, Efficient.h, DefiningPoint.h and Box.h
   * *src* contains DefiningPoint.cpp
  
3) Change in main.cpp line 4 "arguments" to "argument".
    
4) Further instructions can be found in file *build.txt*
