CXX      = $(shell fltk-config --cxx)
DEBUG    = -g
CXXFLAGS = $(shell fltk-config --use-gl --use-images --cxxflags ) -I.
LDFLAGS  = $(shell fltk-config --use-gl --use-images --ldflags )
LDFLAGS += -lstdc++
LDSTATIC = $(shell fltk-config --use-gl --use-images --ldstaticflags )
LINK     = $(CXX)
 
TARGET = fltk_hello
OBJS = fltk_hello.o
SRCS = fltk_hello.cc
 
.SUFFIXES: .o .cc
%.o: %.cc
	$(CXX) $(CXXFLAGS) $(DEBUG) -c $<
 
all: $(TARGET)
	$(LINK) -o $(TARGET) $(OBJS) $(LDSTATIC)
 
$(TARGET): $(OBJS)
fltk_hello.o: fltk_hello.cc
 
clean: $(TARGET) $(OBJS)
	rm -f *.o 2> /dev/null
	rm -f $(TARGET) 2> /dev/null

