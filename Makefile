SRCS=random.cc pri_queue.cc util.cc qalsh.cc srp_lsh.cc l2_alsh.cc \
	l2_alsh2.cc xbox.cc simple_lsh.cc sign_alsh.cc h2_alsh.cc \
	amips.cc pre_recall.cc main.cc
OBJS=${SRCS:.cc=.o}

CXX=g++ -std=c++11
CPPFLAGS=-w -O3

.PHONY: clean

all: ${OBJS}
	${CXX} ${CPPFLAGS} -o alsh ${OBJS}

clean:
	-rm ${OBJS}
