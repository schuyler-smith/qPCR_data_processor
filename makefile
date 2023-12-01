CC:= gcc
CXX:= g++
CXXFLAGS:=
CPPFLAGS= 
LDLIBS= 
LDFLAGS:= -I ./include
HDRS:= $(wildcard include/*) 
SRCS:=

all: \
		bin \
		bin/qPCR_data_processor 

bin:
	$(shell mkdir -p $(@))

bin/qPCR_data_processor: src/qPCR_data_processor.cpp
	$(CXX) $(?) $(CXXFLAGS) -o $(@) $(LDFLAGS)

test: include/test.cpp
	$(CXX) $(?) $(CXXFLAGS) -o $(@) $(LDFLAGS)

# clean:

# install:

.PHONY: bin/qPCR_data_processor

# g++ .\src\qPCR_data_processor.cpp -o qPCR_data_processor -I include -static-libgcc -static-libstdc++
