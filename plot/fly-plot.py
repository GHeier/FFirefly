import plot_path

import argparse
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Process command-line arguments: flags and files."
    )

    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose mode')
    parser.add_argument('--output', type=str, help='Specify output file (use --output=filename)')

    parser.add_argument('files', nargs='*', help='Input file(s)')

    args, unknown = parser.parse_known_args()
    if unknown:
        parser.error(f"Unknown arguments: {' '.join(unknown)}")

    valid_files = []
    for f in args.files:
        if f.startswith('-'):
            parser.error(f"Invalid file name or unknown flag: {f}")
        valid_files.append(f)

    return {
        'flags': {
            'verbose': args.verbose,
            'output': args.output
        },
        'files': valid_files
    }

def get_flags_from_files(result):
    flags = result['flags']
    for filename in result['files']:
        lower_name = filename.lower()
        if 'band' in lower_name:
            flags['band'] = True
        if 'fermi_surface' in lower_name or 'fs' in lower_name or 'FS' in lower_name:
            flags['fermi_surface'] = True

if __name__ == "__main__":
    try:
        result = parse_arguments()
        get_flags_from_files(result)
        #print("Flags:", result['flags'])
        #print("Files:", result['files'])
        if ('band' in result['flags']):
            plot_path.plot_path(result['files'], 'GXMG')

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
