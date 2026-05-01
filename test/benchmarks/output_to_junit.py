import sys
import re
import xml.etree.ElementTree as ET

# This script was AI generated. You have been warned.

def parse_log(log_data):
    total_time = 0.0
    sections = {}

    # 1. Extract total simulation time
    total_match = re.search(r"Simulation finished in ([\d.]+) seconds\.", log_data)
    if total_match:
        total_time = float(total_match.group(1))

    # 2. Extract the Performance Statistics block
    stats_block = ""
    if "=== Performance Statistics ===" in log_data:
        stats_block = log_data.split("=== Performance Statistics ===")[1]
    else:
        # Fallback to searching whole output if the header changes
        stats_block = log_data

    # 3. Dynamically match any section in the format: "  section_name  X.XXX s"
    #    ^\s+      : Leading whitespace
    #    ([\w_]+)  : Capture group 1 (alphanumeric name, e.g., sorting, accumulate_moments)
    #    \s+       : Spaces between name and number
    #    ([\d.]+)  : Capture group 2 (the numeric time, e.g., 12.383)
    #    \s+s$     : Ends with space + 's'
    section_pattern = re.compile(r"^\s+([\w_]+)\s+([\d.]+)\s+s$", re.MULTILINE)
    
    for match in section_pattern.finditer(stats_block):
        # Format key: Replace underscores with spaces and title case it for nice JUnit reports
        name = match.group(1).replace('_', ' ').title()
        time_val = float(match.group(2))
        sections[name] = time_val

    return total_time, sections

def generate_junit_xml(total_time, sections, output_file):
    testsuites = ET.Element("testsuites")
    testsuite = ET.SubElement(testsuites, "testsuite")
    testsuite.set("name", "ThermalCavityBenchmark")
    testsuite.set("tests", str(len(sections) + 1))
    testsuite.set("failures", "0")
    testsuite.set("errors", "0")
    testsuite.set("time", str(total_time))

    # Add overall simulation time as a test case
    tc_total = ET.SubElement(testsuite, "testcase")
    tc_total.set("name", "Total Simulation")
    tc_total.set("classname", "Performance")
    tc_total.set("time", str(total_time))

    # Add dynamically parsed sections
    for name, time_val in sections.items():
        tc = ET.SubElement(testsuite, "testcase")
        tc.set("name", name)
        tc.set("classname", "Performance.Breakdown")
        tc.set("time", str(time_val))

    tree = ET.ElementTree(testsuites)
    tree.write(output_file, encoding="utf-8", xml_declaration=True)

if __name__ == "__main__":
    log_file = sys.argv[1] if len(sys.argv) > 1 else None
    
    if log_file:
        with open(log_file, 'r') as f:
            log_data = f.read()
    else:
        log_data = sys.stdin.read()

    total_time, sections = parse_log(log_data)
    generate_junit_xml(total_time, sections, "benchmark_report.xml")
