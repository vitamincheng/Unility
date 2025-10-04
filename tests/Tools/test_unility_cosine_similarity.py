"""Test functions for utility.py"""

import numpy as np
import pytest
from censo_ext.Tools.utility import cosine_similarity


def test_cosine_similarity_identical_vectors():
    """Test cosine similarity with identical vectors."""
    vec1 = [1, 2, 3]
    vec2 = [1, 2, 3]
    result = cosine_similarity(vec1, vec2)
    assert result == pytest.approx(1.0)


def test_cosine_similarity_opposite_vectors():
    """Test cosine similarity with opposite vectors."""
    vec1 = [1, 2, 3]
    vec2 = [-1, -2, -3]
    result = cosine_similarity(vec1, vec2)
    assert result == pytest.approx(-1.0)


def test_cosine_similarity_orthogonal_vectors():
    """Test cosine similarity with orthogonal vectors."""
    vec1 = [1, 0, 0]
    vec2 = [0, 1, 0]
    result = cosine_similarity(vec1, vec2)
    assert result == pytest.approx(0.0)


def test_cosine_similarity_different_vectors():
    """Test cosine similarity with different vectors."""
    vec1 = [1, 2, 3]
    vec2 = [4, 5, 6]
    result = cosine_similarity(vec1, vec2)
    expected = 0.9746318305482668
    assert result == pytest.approx(expected)


def test_cosine_similarity_zero_vector():
    """Test cosine similarity with zero vector."""
    vec1 = [1, 2, 3]
    vec2 = [0, 0, 0]
    result = cosine_similarity(vec1, vec2)
    assert result == pytest.approx(0.0)


def test_cosine_similarity_different_lengths():
    """Test cosine similarity with different vector lengths."""
    vec1 = [1, 2, 3]
    vec2 = [4, 5]
    with pytest.raises(SystemExit):
        cosine_similarity(vec1, vec2)


def test_cosine_similarity_numpy_arrays():
    """Test cosine similarity with NumPy arrays."""
    vec1 = np.array([1, 2, 3])
    vec2 = np.array([4, 5, 6])
    result = cosine_similarity(vec1, vec2)
    expected = 0.9746318305482668
    assert result == pytest.approx(expected)


def test_cosine_similarity_mixed_types():
    """Test cosine similarity with mixed list and array inputs."""
    vec1 = [1, 2, 3]
    vec2 = np.array([4, 5, 6])
    result = cosine_similarity(vec1, vec2)
    expected = 0.9746318305482668
    assert result == pytest.approx(expected)
