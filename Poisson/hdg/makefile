CFLAGS=-c -Wall -std=c99 -lm -lgsl -lgslcblas
EXEC=ejecutable 

CPPS:=$(wildcard src/*.c)
OBJS:=$(patsubst src/%.c,obj/%.o,$(CPPS))
OBJDIR:=obj

all: bin/$(EXEC)

bin/$(EXEC): $(OBJS) | bin
	gcc -o $@ $(OBJS) -std=c99 -lm -lgsl -lgslcblas
bin: 
	mkdir -p $@

$(OBJDIR)/%.o: src/%.c | $(OBJDIR)
	gcc $(CFLAGS) $< -o $@ 
$(OBJDIR): 
	mkdir -p $@
clean: 
	rm -f -r $(OBJDIR) bin
