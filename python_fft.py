# # import numpy as np
# # import re

# # MAX_N = 256
# # HEADER_REGEX = re.compile(r'^(?P<name>.+\.wav),\s*Min:\s*(?P<min>[-\d]+),\s*Max:\s*(?P<max>[-\d]+)')

# # def apply_hamming_window(data):
# #     """Apply a Hamming window to the data array of length MAX_N."""
# #     return data * np.hamming(MAX_N)

# # def process_block(filename, tokens, stats):
# #     """
# #     Given:
# #       - filename: current .wav block name
# #       - tokens: list of floats [min_placeholder, max_placeholder, sample1, sample2, ...]
# #       - stats: dict tracking global 'max', 'min', and 'count'
# #     Splits samples into 256-point windows, pads the last one if needed,
# #     then for each window applies Hamming + FFT, prints results, and updates stats.
# #     """
# #     # Drop the two placeholders
# #     samples_all = np.array(tokens[2:], dtype=float)
# #     total_samples = samples_all.size

# #     if total_samples == 0:
# #         return

# #     num_windows = (total_samples + MAX_N - 1) // MAX_N

# #     print(f"\nProcessing file: {filename} ({total_samples} samples → {num_windows} windows)")

# #     for w in range(num_windows):
# #         start = w * MAX_N
# #         end = start + MAX_N
# #         window_samples = samples_all[start:end]

# #         # Pad the last window if it's short
# #         if window_samples.size < MAX_N:
# #             window_samples = np.pad(window_samples, (0, MAX_N - window_samples.size), 'constant')

# #         # Apply window and FFT
# #         windowed = apply_hamming_window(window_samples)
# #         fft_res = np.fft.fft(windowed, n=MAX_N)

# #         print(f"\n  Window {w+1}/{num_windows} FFT results:")
# #         for idx, comp in enumerate(fft_res, start=1):
# #             re_val, im_val = comp.real, comp.imag
# #             print(f"    [{idx:3d}] {re_val:.6f} + j{im_val:.6f}")

# #             # Update global min/max
# #             if re_val > stats['max'].real:
# #                 stats['max'] = complex(re_val, stats['max'].imag)
# #             if im_val > stats['max'].imag:
# #                 stats['max'] = complex(stats['max'].real, im_val)
# #             if re_val < stats['min'].real:
# #                 stats['min'] = complex(re_val, stats['min'].imag)
# #             if im_val < stats['min'].imag:
# #                 stats['min'] = complex(stats['min'].real, im_val)
# #             stats['count'] += 1

# # def main():
# #     stats = {
# #         'max': complex(-np.inf, -np.inf),
# #         'min': complex( np.inf,  np.inf),
# #         'count': 0
# #     }

# #     current_name = None
# #     float_tokens = []

# #     with open('min_max_results.txt', 'r') as fp:
# #         for raw_line in fp:
# #             line = raw_line.strip()
# #             if not line:
# #                 continue

# #             # Check for a header line
# #             m = HEADER_REGEX.match(line)
# #             if m:
# #                 # Process the previous file’s block
# #                 if current_name is not None:
# #                     process_block(current_name, float_tokens, stats)

# #                 # Start collecting for the new file
# #                 current_name = m.group('name')
# #                 float_tokens = [
# #                     float(m.group('min')),
# #                     float(m.group('max'))
# #                 ]
# #             else:
# #                 # Not a header: extract floats
# #                 for tok in re.split(r'[\s,]+', line):
# #                     try:
# #                         float_tokens.append(float(tok))
# #                     except ValueError:
# #                         pass

# #         # Process the last file
# #         if current_name is not None:
# #             process_block(current_name, float_tokens, stats)

# #     # Final summary
# #     print("\n=== Summary ===")
# #     mv = stats['max']
# #     nv = stats['min']
# #     print(f"Max number: {mv.real:.6f} + j{mv.imag:.6f}")
# #     print(f"Min number: {nv.real:.6f} + j{nv.imag:.6f}")
# #     print(f"Total FFT outputs: {stats['count']}")

# # if __name__ == "__main__":
# #     main()
# import numpy as np
# import re

# MAX_N = 256

# # This regex captures the header *and* leaves the rest of the line in group 'rest'
# HEADER_REGEX = re.compile(
#     r'^(?P<name>.+?\.wav),\s*Min:\s*(?P<min>-?\d+),\s*Max:\s*(?P<max>-?\d+)\s*(?P<rest>.*)$'
# )

# def apply_hamming_window(data):
#     return data * np.hamming(MAX_N)

# def process_block(filename, tokens, stats):
#     if not tokens:
#         return

#     samples_all = np.array(tokens, dtype=float)
#     total = samples_all.size
#     num_windows = (total + MAX_N - 1) // MAX_N

#     print(f"\n=== {filename} → {total} samples, {num_windows} windows ===")
#     for w in range(num_windows):
#         start = w*MAX_N
#         end = start + MAX_N
#         window = samples_all[start:end]
#         if window.size < MAX_N:
#             window = np.pad(window, (0, MAX_N - window.size), 'constant')

#         windowed = apply_hamming_window(window)
#         fft_res = np.fft.fft(windowed, n=MAX_N)

#         print(f"\n Window {w+1}/{num_windows}:")
#         for i, comp in enumerate(fft_res, 1):
#             r, im = comp.real, comp.imag
#             print(f"  [{i:03d}] {r:.6f} + j{im:.6f}")
#             # update stats
#             if r > stats['max'].real: stats['max'] = complex(r, stats['max'].imag)
#             if im > stats['max'].imag: stats['max'] = complex(stats['max'].real, im)
#             if r < stats['min'].real: stats['min'] = complex(r, stats['min'].imag)
#             if im < stats['min'].imag: stats['min'] = complex(stats['min'].real, im)
#             stats['count'] += 1

# def main():
#     stats = {
#         'max': complex(-np.inf, -np.inf),
#         'min': complex( np.inf,  np.inf),
#         'count': 0
#     }

#     current_name = None
#     buffer = []

#     with open('min_max_results.txt') as f:
#         for raw in f:
#             line = raw.strip()
#             if not line:
#                 continue

#             m = HEADER_REGEX.match(line)
#             if m:
#                 # process previous file
#                 if current_name:
#                     process_block(current_name, buffer, stats)

#                 # start new file
#                 current_name = m.group('name')
#                 # parse min/max placeholders if you want
#                 # ignore them here, so buffer starts fresh
#                 buffer = []

#                 # immediately parse any samples in the 'rest' of the header line
#                 rest = m.group('rest')
#                 for tok in re.split(r'[\s,]+', rest):
#                     try:
#                         buffer.append(float(tok))
#                     except ValueError:
#                         pass
#             else:
#                 # accumulate sample lines
#                 for tok in re.split(r'[\s,]+', line):
#                     try:
#                         buffer.append(float(tok))
#                     except ValueError:
#                         pass

#         # last file
#         if current_name:
#             process_block(current_name, buffer, stats)

#     # summary
#     mv, nv = stats['max'], stats['min']
#     print("\n=== GLOBAL STATS ===")
#     print(f"Max across all: {mv.real:.6f} + j{mv.imag:.6f}")
#     print(f"Min across all: {nv.real:.6f} + j{nv.imag:.6f}")
#     print(f"Total FFT bins: {stats['count']}")

# if __name__ == "__main__":
#     main()
import numpy as np
import re

MAX_N = 256
SAMPLES_PER_FILE = 6000
HEADER_REGEX = re.compile(r'^(?P<name>.+?\.wav),\s*Min:\s*-?\d+,\s*Max:\s*-?\d+')

def apply_hamming_window(data):
    return data * np.hamming(MAX_N)

def process_file(filename, samples, stats):
    total = len(samples)
    num_windows = (total + MAX_N - 1) // MAX_N  # should be 24

    print(f"\n=== {filename}: {total} samples → {num_windows} windows ===")
    for w in range(num_windows):
        start = w * MAX_N
        chunk = samples[start:start+MAX_N]
        if len(chunk) < MAX_N:
            chunk = np.pad(chunk, (0, MAX_N - len(chunk)), 'constant')

        windowed = apply_hamming_window(chunk)
        fft_res = np.fft.fft(windowed, n=MAX_N)

        print(f"\n Window {w+1}/{num_windows}:")
        for i, c in enumerate(fft_res, 1):
            r, im = c.real, c.imag
            print(f"  [{i:03d}] {r:.6f} + j{im:.6f}")
            if r > stats['max'].real: stats['max'] = complex(r, stats['max'].imag)
            if im > stats['max'].imag: stats['max'] = complex(stats['max'].real, im)
            if r < stats['min'].real: stats['min'] = complex(r, stats['min'].imag)
            if im < stats['min'].imag: stats['min'] = complex(stats['min'].real, im)
            stats['count'] += 1

def main():
    stats = {'max': complex(-np.inf, -np.inf),
             'min': complex( np.inf,  np.inf),
             'count': 0}

    with open('min_max_results.txt') as fp:
        while True:
            line = fp.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue

            # Detect header
            m = HEADER_REGEX.match(line)
            if not m:
                continue

            # Found a file header
            filename = m.group('name')
            print(f"\n>>> Found header for {filename}")

            # Now read floats until we have SAMPLES_PER_FILE
            samples = []
            while len(samples) < SAMPLES_PER_FILE:
                data_line = fp.readline()
                if not data_line:
                    break
                for tok in re.split(r'[\s,]+', data_line.strip()):
                    try:
                        samples.append(float(tok))
                        if len(samples) >= SAMPLES_PER_FILE:
                            break
                    except ValueError:
                        pass

            # Process this file
            process_file(filename, np.array(samples, dtype=float), stats)

    # Final summary
    mv, nv = stats['max'], stats['min']
    print("\n=== GLOBAL SUMMARY ===")
    print(f"Max across all bins: {mv.real:.6f} + j{mv.imag:.6f}")
    print(f"Min across all bins: {nv.real:.6f} + j{nv.imag:.6f}")
    print(f"Total FFT bins: {stats['count']}")

if __name__ == "__main__":
    main()
