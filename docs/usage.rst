Usage
=================

Once the PSTools pipeline is correctly initialized, 
you will have two available commands, i.e. pstools and slack.

pstools
------------------

There's an available helper of ``pstools``. One could activate 
it by typing:

.. code-block:: bash

    pstools -h

and return following in terminal

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
                            options:['CI', 'schmidt', 'VST', 'sv', 'REM']
                            (default: None)
      -x XML                xml file that would be served locally (default: None)
      -s {local,eApps,Atlantic_2,Atlantic_3,Linode}
                            GCN server (default: eApps)
      -p {PUBLIC,PRIVATE}   GCN port (default: PUBLIC)
      -v                    Enable task progress report (default: False)
      -l                    Recording task progress report (default: False)

As shown, currenlt there're 5 options from ``pstools``:

- serve alert: 
- Monitor alert:
- Manual input map
- Generate configure file
- Tutorial


slack
--------------------------------

