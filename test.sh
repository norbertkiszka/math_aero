#!/bin/sh

echo "Testing gcc -v
------------------------------------------------";

gcc -v || eval 'echo "gcc non-zero return. Exiting." ; exit 1'

echo "------------------------------------------------";

if [ -f "test" ]; then
	echo "Moving old test executable.";
	echo "mv test test.old";
    mv test test.old;
fi

GCC="gcc -Wall -Werror=format-security -Werror=format -fmax-errors=3 -O2 math_aero.c test.c -o test -lm"

echo $GCC;

$GCC || eval 'echo "Compilation failed. Exiting." ; exit 1'

# echo -e somehow doesnt work properly in my system.

echo "All done. Lets run ./test:
\033[0;32m\033[40m------------------------------------------------\033[0m\e[0m\n";

./test || eval 'echo "\"test\" binary execution failed. Exiting." ; exit 1'

echo "\n\033[0;32m\033[40m------------------------------------------------\033[0m\e[0m\n";

echo "test.sh ended successfully.";

echo "If You found any bug in this library, please report it immediately.";
