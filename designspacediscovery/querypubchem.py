import typing
import requests
import backoff
import re
import time
from tqdm import tqdm
import pickle
import os
import designspacediscovery.utils as ut
import sys
"""
Functions for querying the pubchem pugrest API. 
"""


class pubchemQuery():
    """
    Handles requests to the pubchem PugRest API for use cases in this package. Manage rate limits, timeouts, the like

    No support for asynchronous services
    """

    def __init__(self):
        self.rate_limit_seconds = 5  #requests/second
        self.rate_limit_minutes = 400  # request/minute
        self.request_count_seconds = 0
        self.request_count_minutes = 0
        self.last_second_check = time.time()
        self.last_minute_check = time.time()
        self.count_status = 0
        self.time_status = 0

    def run_queries(self, URLs: dict, cache_params={'cache':False, 'cache_fp':'.', 'cache_name':'cache'}) -> dict:
        """
        For a set of pubchem urls, query pubchem and return the request results. Work on non-batch queries

        Manages pubchem rate limits, retry, other API requesty things for you.

        Parameters:
        ----------
        URLS (dict): dict of {key:URL} pairs
        cache (bool): whether or not to cache intermediate values 
        
        Returns:
        --------

        responses: dict of {key:requests.Response objects or 'FAILED' if issue with request}

        should add something to exit of more than n responses in a row
        """

        assert isinstance(
            URLs, dict), 'This function takes a dictionary of URLS as input'
        assert isinstance(list(URLs.values())[0],
                          str), 'URL in URLs dictionary must be a string'
        

        cache = cache_params['cache']
        cache_fp = cache_params['cache_fp']
        cache_name = cache_params['cache_name']


        response_dict = {}
        print('Querying Pubchem')
        cache_count = 0
        cache_num = 0
        for key in tqdm(list(URLs.keys()), file = sys.stdout, dynamic_ncols=True):
            URL = URLs[key]
            # make sure we are good on pubchem rate limits
            self.__check_rate_status__()
            fail_count = 0
            try:
                response = self.__execute_query__(URL)
                self.__parse_pubchem_header__(response)
            except (requests.exceptions.HTTPError, requests.exceptions.ConnectionError,
         requests.exceptions.ProxyError, requests.exceptions.Timeout,
         requests.exceptions.ReadTimeout) as e:
                tqdm.write(str(e))
                # if the query wrapper has failed, its a lost cause
                response = "FAILED"
                fail_count += 1
            
            response_dict[key] = response
            cache_count +=1

            # quick and dirty caching
            if cache:
                if cache_count > 10:
                    cache_num += 1
                    tqdm.write(f'Pubchem api status: Count status: {self.count_status}, time status: {self.time_status}')
                    with open(f'{cache_fp}/{cache_name}_{cache_num}.pkl', 'wb') as f:
                        pickle.dump(response_dict, f)
                    try:
                        os.remove(f'{cache_fp}/{cache_name}_{cache_num-1}.pkl')
                    except FileNotFoundError as e:
                        pass

                
                cache_count = 0

            # break out if there is something going funky
            if fail_count > 10:
                raise AssertionError('Check URLS or other inputs, 10 failed pubchem requests in a row. Stopping now'
                                     )
        return response_dict
    def batch_queries(self, cid_list, url, max_batch_size=10000):
        """
        Run queries on pubchem endpoints that support batching. This should be most of the pug rest services. Returns a list of the response objects generated by posts, left up to user to decode. 

        Parameters:
        -----------
        cid_list (list of str): list of pubchem cids to query on
        url (str): The url to POST to
        max_batch_size (int) - max size of batches, default 10000

        Returns:
        --------
        response_list: list of response objects

        """
        assert isinstance(
            cid_list, list), 'This function takes a list of cids as input'
        assert isinstance(url, str), 'URL in URLs dictionary must be a string'
        

        print('Querying Pubchem')
        response_list = []
        chunk_num = 0
        for chunk in ut.chunked_iterable(cid_list, max_batch_size):
            #print(chunk)
            print(f'Batch query {100*chunk_num*max_batch_size/len(cid_list):.2f}% complete', end = '\r')
            postbody = 'cid='+','.join(str(cid) for cid in chunk if ut.is_integery(cid))
            self.__check_rate_status__()
            fail_count = 0
            try:
                response = self.__execute_batch_query__(url, postbody)
                self.__parse_pubchem_header__(response)
            except Exception as e:
                print(e)
                # if the query wrapper has failed, its a lost cause
                response = "FAILED"
                fail_count += 1
            
            response_list.append(response)
            chunk_num += 1
            time.sleep(0.5)
        print(f'Batch query 100% complete', end = '\r')
        

                                     
        return response_list

    @backoff.on_exception(
        backoff.expo,
        (requests.exceptions.HTTPError, requests.exceptions.ConnectionError,
         requests.exceptions.ProxyError, requests.exceptions.Timeout,
         requests.exceptions.ReadTimeout),
        max_tries=5)
    def __execute_query__(self, URL: str) -> requests.Response:
        """
        Execute query by URL. Wrapped to enable backoff decorator
        """
        self.request_count_seconds += 1
        self.request_count_minutes += 1
        return requests.get(URL)
    
    @backoff.on_exception(
        backoff.expo,
        (requests.exceptions.HTTPError, requests.exceptions.ConnectionError,
         requests.exceptions.ProxyError, requests.exceptions.Timeout,
         requests.exceptions.ReadTimeout),
        max_tries=5)
    def __execute_batch_query__(self, URL: str, data:dict) -> requests.Response:
        """
        Execute query by URL. Wrapped to enable backoff decorator
        """
        self.request_count_seconds += 1
        self.request_count_minutes += 1
        return requests.post(URL, data = data)

    def __check_rate_status__(self) -> bool:
        """
        check to make sure that we are below the 5 requests/second, and that all the headers check out
        """
        #make sure below requests/second
        if time.time() - self.last_second_check > 5:
            if not self.__second_rate_ok__():
                print('Pubchem requests per second exceeded')
                time.sleep(1)
            else:
                pass
            self.last_second_check = time.time()
        # make sure below requests/minute
        if time.time() - self.last_minute_check > 60:
            if not self.__minute_rate_ok__():
                print(
                    'Pubchem requests per minute exceeded, waiting for 30 seconds'
                )
                time.sleep(30)

        # if count_status or time_status > 75: slight slowdown
        if self.count_status > 75 or self.time_status > 75:
            print('Pubchem status yellow, waiting for a bit')
            time.sleep(10)

        if self.count_status > 95 or self.time_status > 95:
            print('Pubchem requests almost maxed out, waiting')
            time.sleep(30)
        # if count_status or time_status 95: more drastic slowdown

        return

    def __parse_pubchem_header__(self, response):
        header = response.headers['X-Throttling-Control']
        splits = re.split('[()]', header)
        self.count_status = int(splits[1][:-1])
        self.time_status = int(splits[3][:-1])
        return

    def __minute_rate_ok__(self) -> bool:
        """
        Make sure per minute rate limit not exceeded
        """
        rate = (self.request_count_minutes /
                (time.time() - self.last_minute_check)) * 60
        self.request_count_minutes = 0
        return rate < self.rate_limit_minutes

    def __second_rate_ok__(self) -> bool:
        """make sure second rate is ok"""
        rate = self.request_count_seconds / (time.time() -
                                             self.last_second_check)
        self.request_count_seconds = 0
        return rate < self.rate_limit_seconds
