import json

policy_file = "/opt/_internal/pipx/venvs/auditwheel/lib/python3.10/site-packages/auditwheel/policy/manylinux-policy.json"
data = json.load(open(policy_file))
additional_libs = [
        "libbz2.so.1.0.6",
        "libfontconfig.so.1.11.1",
        "libfreetype.so.6.14.0",
        "libGLU.so.1.3.1",
        "libpng15.so.15.13.0",
        "libtcl8.so",
        "libtk8.so",
        "libuuid.so.1.3.0",
        "libz.so.1.2.7",
        "libXmu.so.6",
        "libOpenGL.so.0",
        "libGLdispatch.so.0",
        "libGLX.so.0",
        "libGLU.so.1",
        ]

for entry in data:
    if 'manylinux' in entry['name']:
        entry['lib_whitelist'] += additional_libs

json.dump(data, open(policy_file, 'w'))
