CXX			:= g++
TARGET		:= project.out
BUILDDIR	:= build
SRCDIR		:= src
CXXFLAGS	:= -std=c++17 -g
LIBS		:= -lstdc++
SRCEXT		:= cpp
SOURCES 	:= $(wildcard $(SRCDIR)/*.$(SRCEXT))
OBJECTS		:= $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%, $(SOURCES:.$(SRCEXT)=.o))

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@printf "\e[33m\e[1mBuilding...\e[0m\n";
	@mkdir -p $(BUILDDIR)
	@echo "  $(notdir $@) from $(notdir $<)"
	@$(CXX) $(CXXFLAGS) -c -o $@ $<

$(TARGET): $(OBJECTS)
	@printf "\e[35m\e[1mLinking...\e[0m\n";
	@echo "  $(notdir $(OBJECTS))"
	@$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	@printf "\e[31m\e[1mCleaning...\e[0m\n"
	@echo "  /$(BUILDDIR)"
	@echo "  /$(TARGET)"
	@$(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: r
r: $(TARGET)
	@printf "\e[33m\e[1mRunning $(TARGET)\e[0m\n"
	@./$(TARGET)

.PHONY: run
run: $(TARGET)
	@printf "\e[33m\e[1mRunning $(TARGET)\e[0m\n"
	@./$(TARGET)
