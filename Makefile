SRCS=util.cc pri_queue.cc qalsh_col.cc l2_alsh.cc l2_alsh2.cc xbox.cc \
		h2_alsh.cc srp_lsh.cc simple_lsh.cc sign_alsh.cc amip.cc main.cc
OBJS=$(SRCS:.cc=.o)

CXX=g++ -std=c++11
CPPFLAGS=-w -O3

.PHONY: clean

all: $(OBJS)
	$(CXX) -o alsh $(OBJS)

util.o: util.h

pri_queue.o: pri_queue.h

qalsh_col.o: qalsh_col.h

l2_alsh.o: l2_alsh.h

l2_alsh2.o: l2_alsh2.h

xbox.o: xbox.h

h2_alsh.o: h2_alsh.h

srp_lsh.o: srp_lsh.h

simple_lsh.o: simple_lsh.h

sign_alsh.o: sign_alsh.h

amip.o: amip.h

main.o:

clean:
	-rm $(OBJS) alsh
