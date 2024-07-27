from . import models
from rest_framework import serializers


class SudokuProblemSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = models.Problem
        fields = ['id', 'givens', 'solution', 'difficulty', 'category']
