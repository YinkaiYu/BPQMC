# ==== 自动选择 MPI Fortran 包装器 ====
FC_WRAPPER := $(or \
  $(shell command -v mpiifort 2>/dev/null), \
  $(shell command -v mpiifx   2>/dev/null), \
  $(shell command -v mpifort  2>/dev/null), \
  $(shell command -v mpif90   2>/dev/null), \
  $(shell command -v gfortran 2>/dev/null) \
)
FC_NAME := $(notdir $(FC_WRAPPER))

# 根据包装器选择后端（仅 Intel MPI 有 -fc）
ifeq ($(FC_NAME),mpiifort)
  # 检查ifort是否可用，不可用则使用ifx
  ifeq ($(shell command -v ifort 2>/dev/null),)
    FC_BACKEND := -fc=ifx
  else
    FC_BACKEND := -fc=ifort
  endif
else ifeq ($(FC_NAME),mpiifx)
  FC_BACKEND := -fc=ifx
else
  FC_BACKEND :=
endif

# ==== 设置 HOME，带 fallback ====
HOME_DIR1 := /home/zxli_1/Lib_90_new
HOME_DIR2 := /home/yyk/Lib_90_new

ifeq ($(wildcard $(HOME_DIR1)),)
  HOME := $(HOME_DIR2)
else
  HOME := $(HOME_DIR1)
endif

FC := $(FC_WRAPPER) $(FC_BACKEND)

# ==== 编译/链接参数（标准变量名） ====
FFLAGS := -O3 -c -I$(HOME)/Modules     # 用 -I 而不是复制 .mod
LDFLAGS ?= -mkl
LDLIBS := \
  $(HOME)/Modules/modules_90.a \
  $(HOME)/MyEis/libeis.a \
  $(HOME)/MyNag/libnag.a \
  $(HOME)/MyLin/liblin.a \
  $(HOME)/Ran/libran.a

# 如需保留 SUFFIX 透传
SUFFIX ?=

.PHONY: all clean copymods
.DEFAULT_GOAL := all

# 如确实有代码依赖“当前目录下的 .mod”，保留一个显式目标
copymods:
	@mkdir -p .modcache
	@cp -u $(HOME)/Modules/*.mod .modcache 2>/dev/null || true

all:  ## 调用内层 Compile
	$(MAKE) -f Compile \
	  FC="$(FC)" \
	  LF="$(LDFLAGS)" \
	  FLAGS="$(FFLAGS)" \
	  LIBS="$(LDLIBS)" \
	  SUFFIX="$(SUFFIX)"

clean:
	$(MAKE) -f Compile clean
	@rm -f *.o *.mod *.lst *.opt-report || true
	@rm -rf .modcache || true
