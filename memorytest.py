from partfinder import main, subset, scheme, config

from pympler.classtracker import ClassTracker
tracker = ClassTracker()
# Plug it in here..
scheme.tracker = tracker

scheme.tracker = tracker
tracker.track_class(subset.Subset)
tracker.track_class(scheme.Scheme)
tracker.create_snapshot()
main.call_main("DNA", "--raxml examples/nucleotide --cmd ' -T 2'" )
# main.call_main("DNA", "--raxml Li_2008")
tracker.create_snapshot()
tracker.stats.print_summary()
