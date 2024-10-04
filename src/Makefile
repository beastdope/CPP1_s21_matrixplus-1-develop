CC = g++
CFLAGS = -Wall -Werror -Wextra -c -std=c++17 -pedantic -fsanitize=address -g
CFLAGS_GCOV = -Wall -Werror -Wextra -std=c++17
CHECK_FLAG = -lgtest -lgtest_main -lm -lsubunit 
VALGRIND_FLAGS  = 	--log-file="valgrind.txt" --tool=memcheck --leak-check=yes --track-origins=yes
GCOV_FLAGS = -ftest-coverage -fprofile-arcs -fno-elide-constructors 
TEST_FLAGS = -lgtest -lpthread
TEST = tests/tests.cpp tests/main.cpp 
SRC = $(wildcard *.cpp)
OBJ = $(patsubst %.cpp, %.o, $(SRC))

# -fno-elide-constructors -fno-default-inline to cover the curly brackets {}
SRC_TESTS = $(wildcard tests/*.cpp)
OBJ_TESTS = $(patsubst tests/%.cpp, tests/%.o, $(SRC_TESTS))

all: clean s21_matrix_oop.a

s21_matrix_oop.a: $(OBJ)
	ar cr s21_matrix_oop.a $(OBJ)
test: $(OBJ) $(OBJ_TESTS) 
	$(CC) $(OBJ) $(OBJ_TESTS) $(CHECK_FLAG) -o test
	./test
		
gcov_report : test
	$(CC) $(CFLAGS_GCOV) $(SRC) $(SRC_TESTS) $(CHECK_FLAG) $(GCOV_FLAGS) -o test_coverage
	./test_coverage
	lcov -c -d . --output-file test_coverage.info
	genhtml test_coverage.info --output-directory html_report
	xdg-open ./html_report/index.html
	
valgrind: clean test
	valgrind $(VALGRIND_FLAGS) ./test
	grep -n "ERROR SUMMARY" valgrind.txt

rebuild:
	make clean 
	make all

clean:
	rm -rf *.o
	rm -rf *.a
	rm -rf tests/*.o
	rm -rf *.gcda
	rm -rf *.gcno
	rm -rf *.info
	rm -rf html_report
	rm -rf valgrind.txt
	rm -rf test_coverage
	rm -rf test

clang:
	cp ../materials/linters/.clang-format .
	clang-format -n *.cpp
	clang-format -n */*.cpp
	clang-format -n *.h
	# rm -rf .clang-format

clang-i:
	cp ../materials/linters/.clang-format .
	clang-format -i *.cpp
	clang-format -i */*.cpp
	clang-format -i *.h
	# rm -f .clang-format


