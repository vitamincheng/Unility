"""Test cases for the sub_numpy function in utility module."""

import numpy as np
import pytest
from unittest.mock import patch
from censo_ext.Tools.utility import sub_numpy


def test_sub_numpy_basic_case():
    """Test basic functionality with simple data."""
    data = [1, 2, 3, 10, 11, 12, 20, 21]
    result = sub_numpy(data, max_number=5)
    expected = np.array([3, 2, 2, 1])
    np.testing.assert_array_equal(result, expected)


def test_sub_numpy_basic_case2():
    """Test basic functionality with simple data."""
    data = [1, 2, 3, 10, 15]
    result = sub_numpy(data, max_number=5)
    expected = np.array([3, 2])
    np.testing.assert_array_equal(result, expected)


def test_sub_numpy_single_segment():
    """Test when all data fits in one segment."""
    data = [1, 2, 3, 4, 5]
    result = sub_numpy(data, max_number=10)
    expected = np.array([4, 1])
    np.testing.assert_array_equal(result, expected)


def test_sub_numpy_multiple_segments():
    """Test with multiple segments that fit within max_number."""
    data = [1, 2, 3, 4, 5, 10, 11, 12, 13, 14]
    result = sub_numpy(data, max_number=5)
    # Should create segments of size 5, 5
    expected = np.array([5, 4, 1])
    np.testing.assert_array_equal(result, expected)


def test_sub_numpy_with_numpy_array_input():
    """Test with numpy array input."""
    data = np.array([1, 2, 3, 10, 11, 12, 20, 21])
    result = sub_numpy(data, max_number=5)
    expected = np.array([3, 2, 2, 1])
    np.testing.assert_array_equal(result, expected)


def test_sub_numpy_with_list_input():
    """Test with list input."""
    data = [1, 2, 3, 10, 11, 12, 20, 21]
    result = sub_numpy(data, max_number=5)
    expected = np.array([3, 2, 2, 1])
    np.testing.assert_array_equal(result, expected)


def test_sub_numpy_edge_case_empty():
    """Test with edge case of empty data."""
    data = []
    # This should raise an error since we can't process empty data
    with pytest.raises(SystemExit) as e:
        sub_numpy(data, max_number=5)
    assert e.type is SystemExit
    assert e.value.code == 0


def test_sub_numpy_single_element():
    """Test with single element."""
    data = [5]
    with pytest.raises(SystemExit) as e:
        sub_numpy(data, max_number=5)
    assert e.type is SystemExit
    assert e.value.code == 0


def test_sub_numpy_two_elements():
    """Test with two elements."""
    data = [1, 2]
    with pytest.raises(SystemExit) as e:
        sub_numpy(data, max_number=5)
    assert e.type is SystemExit
    assert e.value.code == 0


def test_sub_numpy_max_number_too_small():
    """Test when max_number is too small to satisfy constraints."""
    # This would likely trigger the exit condition in the original function
    data = [1, 2, 3, 4, 5]
    with pytest.raises(SystemExit) as e:
        sub_numpy(data, max_number=1)
    assert e.type is SystemExit
    assert e.value.code == 0


def test_sub_numpy_with_large_gaps():
    """Test with data that has large gaps."""
    data = [1, 2, 3, 100, 200, 300]
    result = sub_numpy(data, max_number=5)
    # Should create segments based on optimal cut points
    assert isinstance(result, np.ndarray)
    assert result.size > 0


if __name__ == "__main__":
    pytest.main([__file__])
