import typing
import requests
import backoff
import re
import time
"""
Functions for querying the pubchem pugrest API. 
"""

class pubchemQuery():
    """
    Handles requests to the pubchem PugRest API for use cases in this package. Manage rate limits, timeouts, the like

    No support for asynchronous services
    """

    def __init__(self):
        self.rate_limit_seconds = 5 #requests/second
        self.rate_limit_minutes = 400 # request/minute
        self.request_count_seconds = 0
        self.request_count_minutes = 0
        self.last_second_check = time.time()
        self.last_minute_check = time.time()
        self.count_status = 0
        self.time_status = 0


    def run_queries(self, URLs: dict) -> dict:
        """
        For a set of pubchem urls, query pubchem and return the request results.

        Manages pubchem rate limits, retry, other API requesty things for you.

        Parameters:
        ----------
        URLS (dict): dict of {key:URL} pairs
        
        Returns:
        --------

        responses: dict of {key:requests.Response objects}
        """

        assert isinstance(URLs, dict), 'This function takes a dictionary of URLS as input'
        assert isinstance(list(URLs.values())[0], str), 'URL in URLs dictionary must be a string'

        response_dict = {}
        for key in list(URLs.keys()):
            URL = URLs[key]
            # make sure we are good on pubchem rate limits
            self.__check_rate_status__()
            try:
                response = self.__execute_query__(URL)
                self.__parse_pubchem_header__(response)
            except:
                # if the query wrapper has failed, its a lost cause
                response = "FAILED"
            response_dict[key] = response
        return response_dict

    @backoff.on_exception(backoff.expo, (requests.exceptions.HTTPError, requests.exceptions.ConnectionError, requests.exceptions.ProxyError, requests.exceptions.Timeout, requests.exceptions.ReadTimeout), max_tries = 5)
    def __execute_query__(self, URL: str) -> requests.Response:
        """
        Execute query by URL. Wrapped to enable backoff decorator
        """
        self.request_count_seconds += 1
        self.request_count_minutes += 1
        return requests.get(URL)
            
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
                print('Pubchem requests per minute exceeded, waiting for 30 seconds')
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
        rate = self.request_count_minutes/(time.time() - self.last_minute_check)
        self.request_count_minutes = 0
        return rate < self.rate_limit_minutes

    def __second_rate_ok__(self) -> bool:
        """make sure second rate is ok"""
        rate = self.request_count_seconds/(time.time() - self.last_second_check)
        self.request_count_seconds = 0
        return rate < self.rate_limit_seconds