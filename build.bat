clang ^
-D _CRT_SECURE_NO_WARNINGS ^
-D _USE_MATH_DEFINES ^
-std=c99 ^
-O3 ^
-Wall ^
src/*.c ^
-o CFF.exe
CFF.exe examples/ ^
-svg out/ ^
-v out/data.csv
