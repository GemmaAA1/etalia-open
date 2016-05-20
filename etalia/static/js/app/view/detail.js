define([
    'app',
    'text!app/templates/detail.hbs',
    'app/model/detail'
], function (App, template) {

    //var defaults = {};

    App.View.DetailButton = App.Backbone.View.extend({
        tagName: 'div',

        template: App.Handlebars.compile(
            '<button class="btn-circle" type="button">' +
                '<span class="eai eai-{{icon}}"></span>' +
            '</button>' +
            '{{#if caption}}<span>{{caption}}</span>{{/if}}'
        ),

        events: {
            "click": 'onClick'
        },

        render: function() {
            this.$el.html(this.template(this.model.attributes));

            return this;
        },

        onClick: function(e) {
            e.preventDefault();

            this.model.get('callback')(e);
        }
    });

    App.View.Detail = App.Backbone.View.extend({
        tagName: 'div',
        id: 'detail',

        template: App.Handlebars.compile(template),

        events: {
            "click": 'onClick'
        },

        initialize: function (options) {
            if (!options.model) {
                throw 'options.model is mandatory';
            }
        },

        render: function () {
            App.$('div[data-detail-placeholder]').replaceWith(this.$el.html(this.template({})));

            var that = this,
                $bar = this.$('.bar > .wrapper > .inner').empty(),
                $document = this.$('.document > .wrapper').empty();

            App._.forEach(['left', 'right', 'center'], function(name) {
                var button = that.model.get(name + '_button');
                if (button) {
                    var view = new App.View.DetailButton({
                        model: button,
                        attributes: {
                            'class': 'detail-nav detail-nav-' + name,
                            title: button.get('title')
                        }
                    });
                    $bar.append(view.render().$el);
                    that.pushSubView(view);
                }
            });

            var view = this.model.get('view').render();
            $document.append(view.$el);
            that.pushSubView(view);

            this.listenToOnce(view, 'close', function() {
                that.close();
            });

            return this.open();
        },

        onClick: function(e) {
            if (App.$(e.target).attr('id') === 'detail-document') {
                this.close();
            }
        },

        open: function() {
            App.$('body').addClass('detail-opened');

            return this;
        },

        close: function() {
            App.$('body').removeClass('detail-opened');

            this.clearSubViews();

            return this;
        }
    });

    return App.View.Detail;
});
