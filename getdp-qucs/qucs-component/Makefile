PLUGINS=qucsonelab qucsonelab4

CXXFLAGS= -std=gnu++11 -fPIC -shared -rdynamic -DHAVE_CONFIG_H
LDLIBS= -lqucs
PLUGIN_FILES=$(addsuffix .so, $(PLUGINS))

qucsonelab.so: LDLIBS+= -lm
qucsonelab.so: qucsonelab.cpp

qucsonelab4.so: LDLIBS+= -lm
qucsonelab4.so: qucsonelab4.cpp

%.so:
		$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) $^ -o $@

all: $(PLUGIN_FILES)

clean:
		$(RM) $(PLUGIN_FILES)

.PHONY: all clean
.DEFAULT_GOAL=all
