
TIME?=/usr/bin/time

UBMH=$(patsubst bw_HARDLINKED_v8/%,%.bmh.65536.bmh,$(wildcard bw_HARDLINKED_v8/*bw))
LBMH=$(patsubst bw_HARDLINKED/%,%.bmh.65536.bmh,$(wildcard bw_HARDLINKED/*bw))
bmhs=$(UBMH) $(LBMH)
UPMH=$(patsubst bw_HARDLINKED_v8/%,%.pmh.65536.pmh,$(wildcard bw_HARDLINKED_v8/*bw))
LPMH=$(patsubst bw_HARDLINKED/%,%.pmh.65536.pmh,$(wildcard bw_HARDLINKED/*bw))
pmhs=$(UPMH) $(LPMH)
BMH: $(bmhs)
PMH: $(pmhs)
%.bmh.65536.bmh: bw_HARDLINKED_v8/%
	$(TIME) -v ./bwsketch -z65536 -o $(patsubst %.bmh.65536.bmh,%.bmh,$@) $<
%.bmh.65536.bmh: bw_HARDLINKED/%
	$(TIME) -v ./bwsketch -z65536 -o $(patsubst %.bmh.65536.bmh,%.bmh,$@) $<
%.pmh.65536.pmh: bw_HARDLINKED_v8/%
	$(TIME) -v ./bwsketch -z65536 -Po $(patsubst %.pmh.65536.pmh,%.pmh,$@) $<
%.pmh.65536.pmh: bw_HARDLINKED/%
	$(TIME) -v ./bwsketch -z65536 -Po $(patsubst %.pmh.65536.pmh,%.pmh,$@) $<
