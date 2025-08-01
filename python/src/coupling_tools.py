#                                               -*- Python -*-
#
# @brief Gives functions that help coupling against external code,
#   .i.e: manipulate template file.
#
# Copyright 2005-2025 Airbus-EDF-IMACS-ONERA-Phimeca
#
# This library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library.  If not, see <http://www.gnu.org/licenses/>.
#
#

"""
External code helpers.

Provides several functions to ease wrapping of an external code:
- replace: allows one to replace a value in template file
- execute: run an external code
- get: parse values from a result file
"""

import re
import shutil
import os
import shlex
import subprocess
import sys

debug = False
default_encoding = sys.getdefaultencoding()


def check_param(obj, obj_type):
    """Assert obj as type obj_type."""

    try:
        obj_type(obj)
    except Exception:
        raise AssertionError(
            "error: parameter ("
            + str(type(obj))
            + ": "
            + str(obj)
            + ") must be of "
            + str(obj_type)
            + "!"
        )


def replace(infile, outfile, tokens, values, formats=None, encoding=default_encoding):
    """
    Replace values in a file using delimiters.

    Parameters
    ----------
    infile : str
        Name of template file that will be parsed.
    outfile : str
        Name of file that will received the template parsed.
        If equal to *None* or to *infile*, the file will be replaced in-place.
    tokens : list of str
        Regexes that will be replaced.
        Look out for overlapping tokens (eg @x1 and @x10): a solution can be to use padding (@x001, @x010)
        or use markers at both sides (@x1@, @x10@), else longer tokens should be listed first.
    values : list
        Values (can be string, float, ...) that will replace
        the tokens. The list must have the same size as tokens.
    formats : list of str, optional
        A list of formats for printing the values in the parsed files at
        tokens locations
        see https://docs.python.org/2/library/string.html#formatspec
    encoding : str, optional
        File encoding
        see http://docs.python.org/2/library/codecs.html#codec-base-classes

    Raises
    ------
    AssertionError
        parameters badly set
    EOFError
        a token has not been found

    Examples
    --------
    >>> import openturns.coupling_tools as ct
    >>> # write a template file
    >>> with open('prgm.dat.in', 'w') as f:
    ...     count = f.write('E=@E_VAR F=-PF-')
    >>> # replace tokens from template
    >>> ct.replace('prgm.dat.in', 'prgm.dat',
    ...     tokens=["@E_VAR", '-PF-'], values=[1.4, '5'])
    >>> # display file
    >>> with open('prgm.dat', 'r') as f:
    ...     print(f.read())
    E=1.4 F=5
    """
    if len(tokens) != len(values):
        raise AssertionError("Error: tokens size is != values size!")
    check_param(tokens, list)
    check_param(values, list)
    if formats is not None:
        check_param(formats, list)
        if len(values) != len(formats):
            raise AssertionError("Error: values size is != formats size!")
    else:
        formats = ["{0}"] * len(values)

    inplace = False
    if (outfile is None) or (infile == outfile):
        inplace = True
        outfile = infile + ".temporary_outfile"

    regex_tokens = []
    found_tokens = []
    for token in tokens:
        regex_tokens.append(re.compile(token))
        found_tokens.append(False)

    with open(infile, "rb") as fi, open(outfile, "wb") as fo:
        for line in fi:
            line = line.decode(encoding)
            i = 0
            for regex_token in regex_tokens:
                found = regex_token.search(line)
                while found is not None:
                    found_tokens[i] = True
                    line = line.replace(found.group(0), formats[i].format(values[i]))
                    found = regex_token.search(line)
                i += 1
            fo.write(line.encode(encoding))

    for token, found_token in zip(tokens, found_tokens):
        if not found_token:
            raise EOFError(f"No token {token} was found")

    if inplace:
        os.remove(infile)
        shutil.move(outfile, infile)


class OTCalledProcessError(subprocess.CalledProcessError):
    def __str__(self):
        err_msg = (
            (":\n" + self.stderr[:200].decode()) if self.stderr is not None else ""
        )
        return super(OTCalledProcessError, self).__str__() + err_msg


def execute(
    cmd,
    cwd=None,
    shell=False,
    executable=None,
    hide_win=True,
    check=True,
    capture_output=False,
    timeout=None,
    env=None,
):
    """
    Launch an external process.

    Parameters
    ----------
    cmd : list of str or str
        Command to execute, e.g.: "echo 42"
    cwd : str
        Current directory of the executed command.
    shell : bool, default=False
        If set to True, the command is started in a shell (bash).
    executable : str, default=False
        path to the shell. e.g. /bin/zsh.
    hide_win : str, default=True
        Hide cmd.exe popup on windows platform.
    check : bool, default=True
        If set to True: raise RuntimeError if return code of process != 0
    capture_output : bool, default=False
        Whether the output/error streams will be captured
    timeout : int
        Process timeout (Python >=3.3 only)
        On timeout and if psutil is available the children of the process
        are killed before the process itself
    env : dict, default=None
        Environment variables mapping for the new process

    Returns
    -------
    cp : subprocess.CompletedProcess
        Process state info

    Raises
    ------
    RuntimeError
        could not run

    Examples
    --------
    >>> import openturns.coupling_tools as ct
    >>> cp = ct.execute('echo 42', capture_output=True, shell=True)
    >>> cp.returncode
    0
    >>> int(cp.stdout)
    42
    """

    # split cmd if not in a shell before passing it to os.execvp()
    process_args = cmd
    if not shell and isinstance(cmd, str):
        posix = os.name == "posix"
        process_args = shlex.split(cmd, posix=posix)

    # override startupinfo to hide windows console
    startupinfo = None
    if hide_win and sys.platform.startswith("win"):
        startupinfo = subprocess.STARTUPINFO()
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW

    stdout = subprocess.PIPE if capture_output else None
    stderr = subprocess.PIPE if capture_output else None
    process = subprocess.Popen(
        process_args,
        shell=shell,
        cwd=cwd,
        executable=executable,
        stdout=stdout,
        stderr=stderr,
        startupinfo=startupinfo,
        env=env,
    )
    stdout_data = None
    stderr_data = None
    try:
        stdout_data, stderr_data = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        # kill children is psutil is available
        try:
            import psutil
        except ImportError:
            pass
        else:
            parent = psutil.Process(process.pid)
            for child in parent.children(recursive=True):
                child.kill()
        process.kill()
        raise

    returncode = process.poll()

    # check return code
    if check and returncode:
        raise OTCalledProcessError(returncode, cmd, stdout_data, stderr_data)

    return subprocess.CompletedProcess(
        cmd, returncode, stdout=stdout_data, stderr=stderr_data
    )


def get_regex(filename, patterns, encoding=default_encoding):
    """
    Get values from a file using regex.

    Parameters
    ----------
    filename : str
        The name of the file to parse
    patterns : list of str
        Regex patterns that will permit one to get the values
        see https://docs.python.org/2/library/re.html for available patterns
        The value to be searched must be surrounded by parenthesis.
    encoding : str
        File encoding
        see http://docs.python.org/2/library/codecs.html#codec-base-classes

    Returns
    -------
    results : list of float
        Each value corresponds to each pattern.
        If nothing has been found, the corresponding value is set to None.

    Raises
    ------
    AssertionError
        parameters badly set
    EOFError
        no value found

    Examples
    --------
    >>> import openturns.coupling_tools as ct
    >>> # write an output file
    >>> with open('results.out', 'w') as f:
    ...     count = f.write('@E=-9.5E3')
    >>> # parse file with regex
    >>> ct.get_regex('results.out', patterns=[r"(\\R)"])
    [-9500.0]
    """

    if not isinstance(patterns, list) or len(patterns) == 0:
        raise AssertionError("error: patterns parameter badly set!")

    results = [None] * len(patterns)

    re_patterns = []
    for pattern in patterns:
        # OT-like shortcuts
        pattern = pattern.replace(
            r"\R", r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
        )
        pattern = pattern.replace(r"\I", r"[+-]? *\d+")

        re_patterns.append(re.compile(pattern))

    with open(filename, "rb") as f:
        for line in f:
            line = line.decode(encoding)
            i = 0
            for re_pattern in re_patterns:
                match = re_pattern.search(line)
                if match:
                    results[i] = float(match.group(1))
                i += 1

    for result, pattern in zip(results, patterns):
        if result is None:
            raise EOFError(f"No pattern [{pattern}] found")

    return results


def get_real_from_line(line):
    """
    Try to get a real value from the beginning of a text line.

    Parameters
    ----------
    line : str
        Line to parse

    Returns
    -------
    value : float or int
        Found value. Raise an exception if nothing found.
    """

    # \S*: spaces are allowed at the beginning of the real
    real_regex = r"\s*[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"
    re_real = re.compile(real_regex.encode())

    # get the value
    match = re_real.match(line)
    if match:
        result = float(line[match.start() : match.end()])
    else:
        raise EOFError(f"No float found at the beginning of this line: [{line}]")

    return result


def read_line(handle, seek=0):
    """
    Get next line.

    Parameters
    ----------
    seek : int
        if seek < 0: stop reading seek before the end of the file
    """
    line = handle.readline()

    if seek < 0 and len(line) > 0:
        if handle.tell() > -seek:
            if debug:
                sys.stderr.write("line before cut ->" + line.decode() + "<-\n")

            # cut the end of the line
            line = line[: len(line) - (handle.tell() + seek)]
            # after cut, if the line is empty, add a \n to behave like
            # readline() python function
            if len(line) == 0:
                line = os.linesep
            # any further read will return ''
            handle.seek(0, os.SEEK_END)

            if debug:
                sys.stderr.write("line after cut ->" + line.decode() + "<-\n")

    if debug:
        sys.stderr.write("read_line: ->" + line.decode() + "<-\n")

    return line


def get_line_col(
    filename, skip_line=0, skip_col=0, col_sep=None, seek=0, encoding=default_encoding
):
    """
    Get a value at specific line/columns coordinates.

    Parameters
    ----------
    filename : str
        The name of the file that will be parsed
    skip_line : int, default=0
        Number of lines skipped
        If skip_line < 0: count lines backward from the end of the file.
        Be careful: a last empty line is taken into account too!
        Default: 0: no line skipped
    skip_col : int, default=0
        Number of columns skipped from the beginning or end of the line.
        If skip_col < 0: count col backward from the end of the line.
        Default: 0: no column skipped
    col_sep : str
        Column separator
        Default: None: whitespace separator, see str.split
    seek : int, default=0
        if > 0, consider the file starts at pos seek.
        if < 0, consider the file ends at pos -seek (and NOT (end-(-seek))!).
        Default: 0: consider the whole file.
    encoding : str
        File encoding
        see http://docs.python.org/2/library/codecs.html#codec-base-classes

    Returns
    -------
    result : float
        Value to retrieve

    Examples
    --------
    >>> import openturns.coupling_tools as ct
    >>> with open('results.out', 'w') as f:
    ...     count = f.write('1.1 1.2 1.3 1.4')
    >>> ct.get_line_col(filename='results.out', skip_col=2)
    1.3
    """
    check_param(filename, str)
    check_param(skip_line, int)
    check_param(skip_col, int)
    check_param(seek, int)

    with open(filename, "rb") as f:
        if seek > 0:
            f.seek(seek)

        if debug:
            sys.stderr.write(
                "get_line_col(skip_line="
                + str(skip_line)
                + ", skip_col="
                + str(skip_col)
                + ", seek="
                + str(seek)
                + ")\n"
            )

        # skip line backward
        if skip_line < 0:
            # cache position of each beginning of line
            # last elt of the list is the last inserted
            lines_cache = []
            # determine number of elt to cache
            lines_cache_size = -skip_line

            # build lines cache
            previous_pos = f.tell()
            line = read_line(f, seek)
            while len(line) > 0:
                # append to cache
                lines_cache.append(previous_pos)
                if len(lines_cache) > lines_cache_size:
                    # todo: list not optimized to del before the end
                    del lines_cache[0]
                previous_pos = f.tell()

                line = read_line(f, seek)
                if debug:
                    sys.stderr.write("lines_cache: " + str(lines_cache) + "\n")

            if len(lines_cache) < lines_cache_size:
                raise EOFError(f"The file has less than {lines_cache_size} lines")
            else:
                f.seek(lines_cache[0])
                line_found = read_line(f, seek)
                if debug:
                    sys.stderr.write(
                        "line found: ->" + line_found.decode(encoding) + "<-\n"
                    )
        # skip line forward
        else:
            while skip_line >= 0:
                line = read_line(f, seek)
                if skip_line > 0 and len(line) == 0:
                    raise EOFError("The file has less lines than skip_line")
                skip_line -= 1
            line_found = line

    # get the good col
    if skip_col != 0:
        try:
            line_found = line_found.split(
                col_sep.encode() if col_sep is not None else col_sep
            )[skip_col]
        except Exception:
            raise EOFError(
                f"Value not found on this line: [{line_found.decode(encoding)}]"
            )

    # get the value
    result = get_real_from_line(line_found)

    return result


def get_value(
    filename,
    token=None,
    skip_token=0,
    skip_line=0,
    skip_col=0,
    col_sep=None,
    encoding=default_encoding,
):
    """
    Get a value from a file using a delimiter and/or offsets.

    This function is optimized to be rather fast and takes low memory on human
    readable file.

    Parameters
    ----------
    filename : str
        The name of the file that will be parsed
    token : str
        A regex that will be searched.
        The value right after the token is returned.
        Default: None (no token searched)
    skip_token : int
        The number of tokens that will be skipped before getting
        the value. If set to != 0, the corresponding token parameter must
        not be equal to None.
        If skip_tokens < 0: count tokens backward from the end of the file.
        Default: 0: no token skipped
    skip_line : int
        Number of lines skipped from the token found.
        If corresponding token equal None, skip from the beginning of the file.
        If corresponding token != None, skip from the token.
        If skip_line < 0: count lines backward from the token or from the end
        of the file. Be careful: a last empty line is taken into account too.
        Default: 0: no line skipped
    skip_col : int
        Number of columns skipped from the token found.
        If corresponding token = None, skip words from the beginning
        of the line.
        If corresponding token != None, skip words from the token.
        If skip_col < 0: count col backward from the end of the line or from
        the token.
        Default: 0: no column skipped
    col_sep : str
        Column separator
        Default: None: whitespace separator, see str.split
    encoding : str
        File encoding
        see http://docs.python.org/2/library/codecs.html#codec-base-classes

    Returns
    -------
    value : float
        Value found

    Raises
    ------
    AssertionError
        parameters badly set
    EOFError
        no value found

    Examples
    --------
    using a single token

    >>> import openturns.coupling_tools as ct
    >>> with open('results.out', 'w') as f:
    ...     count = f.write('Y1=2.0, Y2=-6.6E56')
    >>> ct.get_value('results.out', token='Y1=')
    2.0

    using token and skip_tokens

    >>> with open('results.out', 'w') as f:
    ...     count = f.write('Y1=2.6 Y1=6.0 # temperature 2')
    >>> ct.get_value('results.out', token='Y1=', skip_token=1)
    6.0

    using column & line

    >>> with open('results.out', 'w') as f:
    ...     count = f.write('1.1 1.2 1.3 1.4')
    >>> ct.get_value(filename='results.out', skip_col=2)
    1.3
    """
    if debug:
        sys.stderr.write(
            "get_value(token="
            + str(token)
            + ", skip_token="
            + str(skip_token)
            + ", skip_line="
            + str(skip_line)
            + ", skip_col="
            + str(skip_col)
            + ")\n"
        )

    # check parameters
    check_param(filename, str)
    if token is not None:
        check_param(token, str)
    elif skip_token != 0:
        raise AssertionError(
            "error: skip_token parameter needs token " "parameter to be set!"
        )
    check_param(skip_token, int)
    check_param(skip_line, int)
    check_param(skip_col, int)

    result = None
    if not token:
        result = get_line_col(
            filename, skip_line, skip_col, col_sep=col_sep, encoding=encoding
        )
    else:
        with open(filename, "rb") as f:

            re_token = re.compile(token.encode(encoding))

            # store previous begin of line pos
            line_pos = f.tell()

            # store token pos [elt 0: start pos of the token, elt 1: end of the
            # token]
            if skip_token >= 0:
                token_pos = None
            else:
                token_pos_cache = []

            line = f.readline()
            while len(line) > 0:
                for token_match in re_token.finditer(line):
                    if skip_token > 0:
                        # token found but it's not the good one
                        skip_token -= 1
                    elif skip_token == 0:
                        # token found
                        token_pos = [
                            line_pos + token_match.start(),
                            line_pos + token_match.end(),
                        ]
                        if debug:
                            sys.stderr.write(
                                "skip_token == 0, line: " + line.decode(encoding) + "\n"
                            )
                        break
                    else:
                        # token wanted in revert order: we first cache them all
                        token_pos_cache.append(
                            [
                                line_pos + token_match.start(),
                                line_pos + token_match.end(),
                            ]
                        )

                if skip_token >= 0 and token_pos is not None:
                    # token found
                    break

                line_pos = f.tell()
                line = f.readline()

            # get the token from the cache
            if skip_token < 0 and len(token_pos_cache) >= -skip_token:
                token_pos = token_pos_cache[skip_token]

            if token_pos is None:
                raise EOFError(f"No token [{token}] found")

            if skip_line == 0 and skip_col == 0:
                # get the real right after the token
                f.seek(token_pos[1])
                line = f.readline()
                if debug:
                    sys.stderr.write(
                        "first token, line_found: " + line.decode(encoding) + "\n"
                    )
                result = get_real_from_line(line)
            else:
                # get the real by skipping from the token
                if skip_line < 0 or (skip_line == 0 and skip_col < 0):
                    # search before the token
                    seek_pos = -token_pos[0]
                    skip_line -= 1  # skip line from the token
                else:
                    # search after the token
                    seek_pos = token_pos[1]
                result = get_line_col(
                    filename,
                    skip_line,
                    skip_col,
                    col_sep=col_sep,
                    seek=seek_pos,
                    encoding=encoding,
                )

    return result


def get(
    filename,
    tokens=None,
    skip_tokens=None,
    skip_lines=None,
    skip_cols=None,
    encoding=default_encoding,
):
    """
    Get several values from a file using delimiters and/or offsets.

    This is equivalent to get_value
    Each following parameters must be a list. The first element of each list
    will be used to get the first value, the 2nd ...
    If the parameter is set to None, the list is set to the default value of
    get_value(), i.e skip_lines=None => [0, 0, ...].

    Parameters
    ----------
    filename : str
        The name of the file that will be parsed
    tokens : list of str
        see get_value
    skip_tokens : list of int
        see get_value
    skip_lines : list of int
        see get_value()
    skip_cols : list of int
        see get_value()
    encoding : str
        File encoding, see
        http://docs.python.org/2/library/codecs.html#codec-base-classes

    Returns
    -------
    result : list of float
        Values parsed

    Raises
    ------
    AssertionError
        parameters badly set
    EOFError
        no value found

    Examples
    --------
    using tokens

    >>> import openturns.coupling_tools as ct
    >>> with open('results.out', 'w') as f:
    ...     count = f.write('Y1=2.0, Y2=-6.6E2')
    >>> ct.get(filename='results.out', tokens=['Y1=', 'Y2='])
    [2.0, -660.0]
    """
    # test parameters and determine the number of value to return
    nb_values = None
    if tokens is not None:
        nb_values = len(tokens)

    err_msg = "error: skip_tokens parameter needs tokens parameter to be set!"
    if skip_tokens is not None:
        if nb_values is None:
            raise AssertionError(err_msg)
        elif nb_values != len(skip_tokens):
            raise AssertionError(err_msg)

    if skip_lines is not None:
        if nb_values is None:
            nb_values = len(skip_lines)
        elif nb_values != len(skip_lines):
            raise AssertionError(err_msg)

    if skip_cols is not None:
        if nb_values is None:
            nb_values = len(skip_cols)
        elif nb_values != len(skip_cols):
            raise AssertionError(err_msg)

    if nb_values is None:
        raise AssertionError("error: no parameters have been set!")

    # init properly every list
    if tokens is None:
        tokens = [None] * nb_values
    if skip_tokens is None:
        skip_tokens = [0] * nb_values
    if skip_lines is None:
        skip_lines = [0] * nb_values
    if skip_cols is None:
        skip_cols = [0] * nb_values

    results = [None] * nb_values

    # possible optimization: get values in one pass
    for n in range(nb_values):
        results[n] = get_value(
            filename,
            tokens[n],
            skip_tokens[n],
            skip_lines[n],
            skip_cols[n],
            encoding=encoding,
        )

    return results
