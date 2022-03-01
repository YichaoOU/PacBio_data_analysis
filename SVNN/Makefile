CC = g++

all: find_SA find_SR combined

find_SA: source/find_SA_in_FastQ.cpp
	${CC} source/find_SA_in_FastQ.cpp -o bin/find_SA_in_FastQ

find_SR: source/find_SR.cpp
	${CC} source/find_SR.cpp -o bin/find_SR

combined: source/combined.cpp
	${CC} source/combined.cpp -o bin/combined