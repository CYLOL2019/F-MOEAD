# Components of Data
This repository contains close prices of stocks and funds.

First, the close prices of all stocks and funds from 20091211 to 20201211 are downloaded and preprocessed.

Second, we set the time interval as 730 days and the time window as 365 days to reorganize the data for our experiments.
Hence, data in 20091211-20111211, 20101211-20121210,..., 20181209-20201208 are documented.

Finally, according to our needs, 10 instances are picked. They are listed as follows
Time                            #Assets             Component(s)
20091211-20111211   475                       funds
20101211-20121210   595                       funds
20111211-20131210   776                       funds
20091211-20111211   1601                     stocks
20171209-20191209   3406                     stocks
20181209-20201208   3515                     stocks
20171209-20191209   4558                     funds
20181209-20201208   5575                     funds
20171209-20191209   7964                     stocks+funds
20181209-20201208   9090                     stocks+funds

Only the data for the first test problem is uploaded due to the space limitation of Github.
