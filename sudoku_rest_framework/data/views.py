from rest_framework import permissions, viewsets

from . import models
from . import serializers


class SudokuViewSet(viewsets.ModelViewSet):
    """
    API endpoint that allows users to be viewed or edited.
    """
    queryset = models.Problem.objects.all()
    serializer_class = serializers.SudokuProblemSerializer
