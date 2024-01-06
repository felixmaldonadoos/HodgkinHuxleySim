import subprocess

# Path to the compiled C++ executable
executable_path = "./example"  # Change this to the path of your executable

# Run the executable and capture its output
try:
    result = subprocess.run(executable_path, capture_output=True, text=True, check=True)
    print("Output from C++ program:")
    print(result.stdout)
except subprocess.CalledProcessError as e:
    print("Error running C++ program:", e)
