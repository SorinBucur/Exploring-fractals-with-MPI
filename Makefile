all: build

build: main.cpp
	mpic++ main.cpp -o main

clean:
	rm main
