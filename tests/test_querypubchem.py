import designspacediscovery.querypubchem as qpc
import requests
from unittest.mock import patch
import time


class Test_pubchemQuery:

    def test_parse_pubchem_header(self):
        pcquery = qpc.pubchemQuery()
        response = requests.Response()
        response.headers = {'X-Throttling-Control':'Request Count status: Green (0%), Request Time status: Green (0%), Service status: Green (10%)'}
        # test that the pubchem header parser correctly pulls out rate stuff
        pcquery.__parse_pubchem_header__(response)

        assert pcquery.count_status == 0, 'Error in extracting count status from header'
        assert pcquery.time_status == 0, 'Error in extracting time status'

    @patch('designspacediscovery.querypubchem.requests.get')
    def test__execute_query(self, mock_get):
        pcquery = qpc.pubchemQuery()
        url = 'abcd'
        mock_get.return_value.ok = True

        response = pcquery.__execute_query__(url)

        assert pcquery.request_count_minutes == 1
        assert pcquery.request_count_seconds == 1
        
        # 1. Test that request counts are incremented correctly
    @patch('designspacediscovery.querypubchem.requests.get')
    def test_run_queries(self, mock_get):
        pcquery = qpc.pubchemQuery()
        url = 'abcd'
        #1. test that if get function mocked to return an HTTPError, there are 5 requests
        mock_get.side_effect = requests.exceptions.HTTPError
        pcquery.run_queries({'q':url})
        assert pcquery.request_count_seconds == 5

    def test__check_rate_status(mock_sleep):
        # this implicitly tests the rate checker functions (min/sec_rate_ok__()) as well
        pcquery = qpc.pubchemQuery()
        pcquery.last_minute_check = time.time() - 61
        pcquery.last_second_check = time.time()

        # if second rate too high, waits 1 second
        #use patch context so runner doesn't have to wait for big sleep
        with patch('designspacediscovery.querypubchem.time.sleep', return_value = None) as mock_time:
            pcquery.request_count_minutes = 500
            pcquery.__check_rate_status__()
            assert mock_time.called

        # if minute rate is too high, waits 30 seconds
        pcquery.last_minute_check = time.time() 
        pcquery.last_second_check = time.time() - 5.5
        with patch('designspacediscovery.querypubchem.time.sleep', return_value = None) as mock_time:
            pcquery.request_count_seconds = 35 # if we've made 35 requests in the last ~6 seconds we're over limit
            pcquery.__check_rate_status__()
            assert mock_time.called

        # if count statuses > 75, waits 10 seconds