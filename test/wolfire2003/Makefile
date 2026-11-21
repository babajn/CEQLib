CC ?= cc
GSL_DIR ?= $(or $(shell gsl-config --prefix 2>/dev/null),$(wildcard $(HOME)/lib/gsl),/usr)
EQCLIB_DIR := ../../src
OUT_DIR := ../out/wolfire2003
WRAPPED_RATES := SoftXrayIonizationRateHI SoftXrayIonizationRateHeI \
                 TotalSoftXrayHeatingRate CosmicRayIonizationRateHI TotalCosmicRayHeatingRate
comma := ,

CPPFLAGS += -I$(EQCLIB_DIR) -I$(GSL_DIR)/include
CFLAGS += -std=c99 -O3 -Wall -Wextra
LDFLAGS += $(addprefix -Wl$(comma)--wrap=,$(WRAPPED_RATES)) -Wl,-rpath,$(GSL_DIR)/lib
LDLIBS += -L$(GSL_DIR)/lib -lgsl -lgslcblas -lm

.PHONY: run eqclib
run: $(OUT_DIR)/wolfire2003.out
	@set -e; csv=$$(mktemp "$(OUT_DIR)/equilibrium.XXXXXX"); \
	trap 'rm -f "$$csv"' EXIT; \
	$(OUT_DIR)/wolfire2003.out > "$$csv"; \
	mv "$$csv" "$(OUT_DIR)/equilibrium.csv"

eqclib:
	$(MAKE) -C $(EQCLIB_DIR) GSL_DIR=$(GSL_DIR) libEQCLib.a

$(OUT_DIR)/wolfire2003.out: Wolfire2003.c eqclib
	mkdir -p $(OUT_DIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $< $(EQCLIB_DIR)/libEQCLib.a $(LDLIBS) -o $@
