Usage
=================

Once the PSTools pipeline is correctly initialized, 
you will have two available commands, i.e. pstools 
and slack, see `procedures <https://pstools-documentation.readthedocs.io/en/latest/ref.html>`_.

pstools
------------------

There's an available helper of ``pstools``. One could activate 
it by typing:

.. code-block:: bash

    pstools -h

and would return the following lines

.. code-block:: bash

    usage: pstools.py [-h] [-f FITS] [-t TEL] [-x XML]
                  [-s {local,eApps,Atlantic_2,Atlantic_3,Linode}]
                  [-p {PUBLIC,PRIVATE}] [-v] [-l]
                  {1,2,3,4,5}

    >> PStools main algorithm

    positional arguments:
      {1,2,3,4,5}           PSTools mode: [1] Serve alert; [2] Monitor alert; [3]
                            Manual input map; [4] Generate configure file; [5]
                            Tutorial

    optional arguments:
      -h, --help            show this help message and exit
      -f FITS               handle healpix fits file (default: None)
      -t TEL                specific telescopes, if multiple divided with `,`
                            options:['VST', 'REM', ... ...]
                            (default: None)
      -x XML                xml file that would be served locally (default: None)
      -s {local,eApps,Atlantic_2,Atlantic_3,Linode}
                            GCN server (default: eApps)
      -p {PUBLIC,PRIVATE}   GCN port (default: PUBLIC)
      -v                    Enable task progress report (default: False)
      -l                    Recording task progress report (default: False)

As shown, there're 5 options for ``pstools``:

- Serve alert: For a test purpose, you can serve xml via the local server, afterwards hire PSTools to monitor the local server alerts.

- Monitor alert: You can activate PSTools to monitor a server/port, so that, the sudden alert would activate the prioritization process.

- Manual input map: This option would ask users to input a trigger map in healpix format, and the prioritization was subsequently performed.

- Generate configure file: Since there're so many parameters and options in the process, in order to provide users a clear view, let users fast revise them, there're a list of configure files, which stored the parameters, would be used. There's a general configure file, named pstools.default, together with few telescope configure files, named pst_[`tel name`].default.

- Tutorial: will show the tutrial.


slack
--------------------------------

to be done
