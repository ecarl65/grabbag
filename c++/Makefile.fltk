CXX      = $(shell fltk-config --cxx)
DEBUG    = -g
CXXFLAGS = $(shell fltk-config --use-gl --use-images --cxxflags ) -I.
LDFLAGS  = $(shell fltk-config --use-gl --use-images --ldflags ) -lstdc++
LDSTATIC = $(shell fltk-config --use-gl --use-images --ldstaticflags )
LINK     = $(CXX)
 

all: fltk_hello Button

fltk_hello: fltk_hello.o
	$(CXX) $(CXXFLAGS) $(DEBUG) $< -o $@ $(LDFLAGS)

Button: Button.o
	$(CXX) $(CXXFLAGS) $(DEBUG) $< -o $@ $(LDFLAGS)

fltk_hello.o: fltk_hello.cc
	$(CXX) $(CXXFLAGS) $(DEBUG) -c $<

Button.o: Button.cc
	$(CXX) $(CXXFLAGS) $(DEBUG) -c $<
 
clean: $(TARGET) $(OBJS)
	rm -f *.o 2> /dev/null
	rm -f $(TARGET) 2> /dev/null

