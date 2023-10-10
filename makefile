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
		bin/smart_chip_analyzer 

bin:
	$(shell mkdir -p $(@))

bin/smart_chip_analyzer: src/smart_chip_analyzer.cpp
	$(CXX) $(?) $(CXXFLAGS) -o $(@) $(LDFLAGS)

test: include/test.cpp
	$(CXX) $(?) $(CXXFLAGS) -o $(@) $(LDFLAGS)

# clean:

# install:

.PHONY: bin/smart_chip_analyzer

# g++ .\src\smart_chip_analyzer.cpp -o smartchip_analyzer_win -I include -static-libgcc -static-libstdc++
