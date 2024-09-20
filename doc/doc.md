## Documentation Editing

For help editing the documentation visit [mkdocs.org](https://www.mkdocs.org). To generate the docs locally type in the parent directory: `mkdocs serve`
and point the browser to [127.0.0.1.8000](http://127.0.0.1:8000)

You will need to install the `python-markdown-math` extension for rendering equations and the `markdown-callouts` extension for correctly displaying the warning and note blocks. All requirements can be installed automatically using

```bash
$ pip install -r docs/requirements.txt
```
To  install **mkdocs** 

```bash
$ pip install mkdocs
```
or dependencies `$ pip install mkdocs`



```bash
$ pip install pip-tools
```

if you add new markdown extensions, edit the `requirements.in`  file under `docs/`

```bash
$ pip-compile requirements.in
```




