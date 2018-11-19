{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block methods %}
      .. automethod:: __init__

      {% if methods %}
         .. rubric:: Methods

         .. autosummary::
            :toctree:
            {% for item in methods %}
               {%- if not item.startswith('_') %}
                  ~{{ name }}.{{ item }}
               {%- endif -%}
            {%- endfor %}
      {% endif %}
   {% endblock %}

   {% block attributes %}
      {% if attributes %}
         .. rubric:: Attributes

         .. autosummary::
            :toctree:
            {% for item in attributes %}
               ~{{ name }}.{{ item }}
            {%- endfor %}
      {% endif %}
   {% endblock %}
