CXX ?= g++
CXXFLAGS ?= -O3 -DNDEBUG -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS += -Iinclude
LDLIBS += -lntl -lgmp -pthread
ifneq ($(NTL_PREFIX),)
CPPFLAGS += -I$(NTL_PREFIX)/include
LDFLAGS += -L$(NTL_PREFIX)/lib -Wl,-rpath,$(NTL_PREFIX)/lib
endif

.PHONY: all clean
all: bin/hhy_bench

bin/hhy_bench: src/main.cpp $(wildcard include/*.hpp)
	mkdir -p bin
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/main.cpp $(LDFLAGS) $(LDLIBS) -o $@

clean:
	rm -f bin/hhy_bench
