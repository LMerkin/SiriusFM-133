#If some problems with Makefile are occured, consider running with
#g++ -Wall -g -std=c++17 Test2.cpp IRProvider.cpp -o a && ./a [args]

TARGET  = ../Test5
SOURCES = Test5 IRProviderConst 

CXX     = g++
CXXFLAGS += -fopenmp
EXTLIBS = -lgomp

#CXX      = nvc++
#CXXFLAGS += -acc=gpu # -Minfo messages
#EXTLIBS = -lacchost

CXXFLAGS  += -std=c++17 -Wall -Wno-stringop-truncation
CXXFLAGS  += -O3 -DNDEBUG -march=native
#CXXFLAGS += -O0 -g

BUILD_DIR = $(shell pwd)
OBJECTS_DIR = $(BUILD_DIR)/../obj
OBJECTS = $(patsubst %, $(OBJECTS_DIR)/%.o, $(SOURCES))

all: $(TARGET)

$(OBJECTS_DIR):
	$(shell mkdir -p $(OBJECTS_DIR))

$(OBJECTS_DIR)/%.o: %.cpp 
	@echo [CXX] $(@F)
	$(CXX) $(CXXFLAGS) -o $@ -c $(realpath $<)

$(OBJECTS_DIR)/%.o: %.c
	@echo [CC] $(@F)
	$(CC) $(CFLAGS) -o $@ -c $(realpath $<)

$(TARGET) : $(OBJECTS_DIR) $(OBJECTS)
	$(CXX) -o $(TARGET) $(OBJECTS) $(EXTLIBS)

clean:
	$(shell rm -fr $(OBJECTS_DIR))
	$(shell rm -f $(TARGET))
