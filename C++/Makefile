TOP:= $(realpath ./)
SOURCES:=${wildcard $(TOP)/*.cpp}
SOURCES+=${wildcard $(TOP)/geometries/*.cpp}

BUILDDIR:=build

OBJECTS:= $(addprefix $(TOP)/$(BUILDDIR)/,$(subst $(TOP)/,,$(SOURCES:%.cpp=%.o)))

main: $(OBJECTS)
	$(CXX) -std=c++17 -pthread -o $@ $^

$(OBJECTS): $(TOP)/$(BUILDDIR)/%.o: $(TOP)/%.cpp
	$(CXX) -std=c++17 -pthread -c $< -o $@

clean:
	@rm -f $(TOP)/$(BUILDDIR)/*.o