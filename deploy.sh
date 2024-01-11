./lint.sh && \
./test.sh && \
python3 -m pip install --upgrade build && \
python3 -m pip install --upgrade twine && \
rm -rf dist && \
python3 -m build && \
python3 -m twine upload dist/*
