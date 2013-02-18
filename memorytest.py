from partfinder import main, subset, scheme, config

from pympler.classtracker import ClassTracker
tracker = ClassTracker()
scheme.tracker = tracker
tracker.track_class(subset.Subset)
tracker.track_class(scheme.Scheme)
tracker.create_snapshot()
main.call_main("DNA", "--force --raxml --cmdline-extras ' -T 2 ' examples/nucleotide")
tracker.create_snapshot()
main.call_main("DNA", "--raxml --cmdline-extras ' -T 2 ' examples/nucleotide")
tracker.create_snapshot()
tracker.stats.print_summary()
