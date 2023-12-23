TOP:= $(realpath ./)
SOURCES:=${wildcard $(TOP)/main.cpp}
SOURCES+=${wildcard $(TOP)/geometries/*.cpp}
SOURCES := $(filter-out $(TOP)/pybind.cpp, $(SOURCES))

BUILDDIR:=build

OBJECTS:= $(addprefix $(TOP)/$(BUILDDIR)/,$(subst $(TOP)/,,$(SOURCES:%.cpp=%.o)))

main: $(OBJECTS)
	$(CXX) -std=c++17 -pthread -g -o $@ $^

$(OBJECTS): $(TOP)/$(BUILDDIR)/%.o: $(TOP)/%.cpp
	$(CXX) -std=c++17 -pthread -g -c $< -o $@

clean:
	@rm -f $(TOP)/$(BUILDDIR)/*.o