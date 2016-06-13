define([
    'app',
    'text!app/templates/ui/controls.hbs',
    'text!app/templates/ui/controls/search.hbs',
    'text!app/templates/ui/controls/cluster.hbs',
    'text!app/templates/ui/controls/timespan.hbs'
], function (App, template, searchTemplate, clusterTemplate, timespanTemplate) {

    var ControlModel = App.Backbone.Model.extend({
        defaults: {
            active: true,
            value: null
        }
    });

    var ControlView = App.Backbone.View.extend({
        controlsView: null,

        initialize: function(options) {
            if (!options.controlsView) {
                throw 'options.controlsView is mandatory';
            }
            this.controlsView = options.controlsView;

            if (!options.model) {
                this.model = new ControlModel();
            }

            this.listenTo(this.model, 'change', this.onChange);
        },

        onChange: function() {
            this._handleVisibility();

            this.controlsView.dispatchChange();
        },

        _handleVisibility: function() {
            if (this.model.get('active')) {
                this._show();
            } else {
                this._hide();
            }
        },

        _show: function() {
            this.$el.show();
        },

        _hide: function() {
            this.$el.hide();
        },

        enable: function() {
            this.model.set('active', true);
        },

        disable: function() {
            this.model.set('active', false);
        },

        render: function () {
            this.$el.html(this.template({}));

            this._handleVisibility();

            return this;
        }
    });

    var SelectorView = ControlView.extend({
        events: {
            "click button": "onButtonClick",
            "click .choices a": "onChoiceClick"
        },

        open: function() {
            this.$el.addClass('opened');
        },

        close: function() {
            this.$el.removeClass('opened');
        },

        onButtonClick: function(e) {
            e.preventDefault();

            if (this.$el.hasClass('opened')) {
                this.close();
            } else {
                this.open();
            }
        },

        onChoiceClick: function(e) {
            e.preventDefault();

            this.model.set('value', App.$(e.target).closest('a').data('value'));
        }
    });

    var ClusterView = SelectorView.extend({
        tagName: 'div',
        attributes: {
            id: 'cluster',
            'class': 'selector',
            title: 'Clusters'
        },
        template: App.Handlebars.compile(clusterTemplate),

        values: {
            0: {label: 'All'},
            1: {label: 'Purple'},
            2: {label: 'Red'},
            3: {label: 'Orange'},
            4: {label: 'White'}
        },

        _toggleIconVisibility: function(value) {
            if (0 < value) {
                this.$('#cluster-selection').css({display: 'inline-block'});
                this.$('#cluster-none').hide();
            } else {
                this.$('#cluster-selection').hide();
                this.$('#cluster-none').css({display: 'inline-block'});
            }
        },

        onChange: function() {
            this._handleVisibility();

            var value = this.model.get('value');
            if (App._.has(this.values, value)) {
                this.$('#cluster-selection').removeClass('cluster-0 cluster-1 cluster-2 cluster-3');
                if (0 < value) {
                    this.$('#cluster-selection').addClass('cluster-' + (value - 1));
                }
                this._toggleIconVisibility(value);
                this.close();
                this.controlsView.dispatchChange();
            }
        },

        render: function () {
            var value = this.model.get('value');

            this.$el.html(this.template({
                cluster: value - 1
            }));

            this._toggleIconVisibility(value);

            this._handleVisibility();

            return this;
        }
    });

    var TimespanView = SelectorView.extend({
        tagName: 'div',
        attributes: {
            id: 'timespan',
            'class': 'selector',
            title: 'Time span'
        },
        template: App.Handlebars.compile(timespanTemplate),

        values: {
            7:   {icon: 'W',  label: 'Week',     html: '<strong>W</strong>eek'},
            30:  {icon: '1m', label: '1 month',  html: '<strong>1 m</strong>onth'},
            60:  {icon: '2m', label: '2 months', html: '<strong>2 m</strong>onths'},
            365: {icon: '1y', label: '1 year',   html: '<strong>1 y</strong>ear'}
        },

        onChange: function() {
            this._handleVisibility();

            var value = this.model.get('value');
            if (App._.has(this.values, value)) {
                this.$('#timespan-selection').html(this.values[value].icon);
                this.close();
                this.controlsView.dispatchChange();
            }
        },

        render: function () {
            this.$el.html(this.template({
                value: this.values[this.model.get('value')].icon,
                values: this.values
            }));

            this._handleVisibility();

            return this;
        }
    });

    var SearchView = ControlView.extend({
        tagName: 'div',
        attributes: {
            id: 'search'
        },
        template: App.Handlebars.compile(searchTemplate),

        events: {
            "keyup input": "onInputKeyUp",
            "click button[type=submit]": "onSubmitClick",
            "click #close-search": "onCloseClick"
        },

        keyUpTimeout: null,

        _handleVisibility: function() {},

        _clearTimeout: function() {
            if (this.keyUpTimeout) {
                clearTimeout(this.keyUpTimeout);
            }
        },

        open: function() {
            this.$el.addClass('opened');
        },

        close: function() {
            this.$el.removeClass('opened');
        },

        onInputKeyUp: function() {
            this._clearTimeout();

            var model = this.model,
                $input = this.$('input');
            this.keyUpTimeout = setTimeout(function() {
                model.set('value', $input.val())
            }, 800);
        },

        onSubmitClick: function() {
            this._clearTimeout();

            this.model.set('value', this.$('input').val())
        },

        onCloseClick: function() {
            this.close();
        }
    });

    var PinView = ControlView.extend({
        tagName: 'button',
        attributes: {
            id: 'toggle-pinned',
            'class': 'btn-circle',
            type: 'button',
            title: 'Pinned items'
        },

        events: {
            "click": "onClick"
        },

        onClick: function(e) {
            e.preventDefault();

            this.model.set('value', !this.model.get('value'));
        },

        onChange: function() {
            this._handleVisibility();

            if (this.model.get('value')) {
                this.$el.addClass('active');
            } else {
                this.$el.removeClass('active');
            }

            this.controlsView.dispatchChange();
        },

        render: function() {
            this.$el.html('<span class="eai eai-pin"></span>');

            this._handleVisibility();

            return this;
        }
    });

    App.View.Ui = App.View.Ui || {};

    App.View.Ui.Controls = App.Backbone.View.extend({
        tagName: 'div',
        className: 'wrap',

        template: App.Handlebars.compile(template),

        search: null,
        cluster: null,
        timespan: null,
        pin: null,

        events: {
            "click #toggle-search": "onToggleSearchClick"
        },

        initialize: function() {
            this.search = new SearchView({
                controlsView: this
            });
            this.cluster = new ClusterView({
                controlsView: this
            });
            this.timespan = new TimespanView({
                controlsView: this,
                model: new ControlModel({
                    value: 7
                })
            });
            this.pin = new PinView({
                controlsView: this
            });
        },

        onToggleSearchClick: function() {
            this.search.open();
        },

        getContext: function() {
            var data = {};

            if (true === this.pin.model.get('value')) {
                data.pinned = 1;
            }

            // TODO Timespan

            // TODO Cluster

            var search = this.search.model.get('value');
            if (search && 0 < String(search).length) {
                data.search = search;
            }

            return data;
        },

        dispatchChange: function() {
            this.trigger('context-change');
        },

        render: function () {
            this.$el.html(this.template({}));

            this.$('[data-search-placeholder]').replaceWith(this.search.render().$el);
            this.pushSubView(this.search);

            this.$('[data-cluster-placeholder]').replaceWith(this.cluster.render().$el);
            this.pushSubView(this.cluster);

            this.$('[data-timespan-placeholder]').replaceWith(this.timespan.render().$el);
            this.pushSubView(this.timespan);

            this.$('[data-pin-placeholder]').replaceWith(this.pin.render().$el);
            this.pushSubView(this.pin);

            return this;
        }
    });

    App.View.Ui.Controls.create = function(options, createOptions) {
        options = options || {};

        var view = new App.View.Ui.Controls(options);
        if (createOptions) {
            App.View.create(view, createOptions);
        }

        return view;
    };

    return App.View.Ui.Controls;
});
