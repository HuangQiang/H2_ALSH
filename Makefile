SRCS=random.cc pri_queue.cc util.cc qalsh.cc srp_lsh.cc l2_alsh.cc \
	l2_alsh2.cc xbox.cc simple_lsh.cc sign_alsh.cc h2_alsh.cc \
	amips.cc pre_recall.cc main.cc
OBJS=${SRCS:.cc=.o}

CXX=g++ -std=c++11
CPPFLAGS=-w -O3

.PHONY: clean

all: ${OBJS}
	${CXX} ${CPPFLAGS} -o alsh ${OBJS}

random.o: random.h

pri_queue.o: pri_queue.h

util.o: util.h

qalsh.o: qalsh.h

srp_lsh.o: srp_lsh.h

l2_alsh.o: l2_alsh.h

l2_alsh2.o: l2_alsh2.h

xbox.o: xbox.h

simple_lsh.o: simple_lsh.h

sign_alsh.o: sign_alsh.h

h2_alsh.o: h2_alsh.h

amips.o: amips.h

pre_recall.o: pre_recall.h 

main.o:

clean:
	-rm ${OBJS}
